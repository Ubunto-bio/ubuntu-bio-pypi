"""Offline Natural Earth 1:50m polygon validation for African Union member states."""

from __future__ import annotations

import functools
import gzip
import json
import math
from pathlib import Path
from typing import Final

DEFAULT_COASTAL_BUFFER_DEG: Final[float] = 0.25
_BOUNDARIES_PATH: Final[Path] = Path(__file__).with_name("au_boundaries.json.gz")

PolygonRings = list[list[list[float]]]
CountryBBox = tuple[float, float, float, float]


@functools.lru_cache(maxsize=1)
def _load_country_polygons() -> dict[str, tuple[list[PolygonRings], CountryBBox]]:
    """Load quantized Natural Earth 1:50m AU boundaries and precompute bounding boxes."""
    raw: dict[str, list[PolygonRings]] = json.loads(
        gzip.decompress(_BOUNDARIES_PATH.read_bytes()).decode("utf-8")
    )
    indexed: dict[str, tuple[list[PolygonRings], CountryBBox]] = {}
    for code, polys in raw.items():
        min_lon = 180.0
        min_lat = 90.0
        max_lon = -180.0
        max_lat = -90.0
        for poly in polys:
            for ring in poly:
                for lon, lat in ring:
                    if lon < min_lon:
                        min_lon = lon
                    if lon > max_lon:
                        max_lon = lon
                    if lat < min_lat:
                        min_lat = lat
                    if lat > max_lat:
                        max_lat = lat
        indexed[code] = (polys, (min_lon, min_lat, max_lon, max_lat))
    return indexed


def _point_in_ring(lon: float, lat: float, ring: list[list[float]]) -> bool:
    """Return True when point lies inside ring using even-odd ray casting."""
    inside = False
    n = len(ring)
    if n < 3:
        return False
    j = n - 1
    for i in range(n):
        xi, yi = ring[i]
        xj, yj = ring[j]
        if (yi > lat) != (yj > lat):
            denom = yj - yi
            x_intersect = (xj - xi) * (lat - yi) / denom + xi
            if lon < x_intersect:
                inside = not inside
        j = i
    return inside


def _point_in_polygon(lon: float, lat: float, poly: PolygonRings) -> bool:
    """Return True when point lies inside exterior ring and outside interior holes."""
    if not poly or not _point_in_ring(lon, lat, poly[0]):
        return False
    return all(not _point_in_ring(lon, lat, hole) for hole in poly[1:])


def _point_to_segment_dist_sq(
    px: float,
    py: float,
    ax: float,
    ay: float,
    bx: float,
    by: float,
) -> float:
    """Compute squared Euclidean distance in degrees from point (px, py) to segment AB."""
    dx = bx - ax
    dy = by - ay
    len_sq = dx * dx + dy * dy
    if len_sq == 0.0:
        vx, vy = px - ax, py - ay
        return vx * vx + vy * vy
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / len_sq))
    proj_x = ax + t * dx
    proj_y = ay + t * dy
    rx = px - proj_x
    ry = py - proj_y
    return rx * rx + ry * ry


def _within_coastal_buffer(
    lon: float,
    lat: float,
    polys: list[PolygonRings],
    buffer_deg: float,
) -> bool:
    """Return True when point is within buffer_deg of any exterior boundary segment."""
    if buffer_deg <= 0.0:
        return False
    max_sq = buffer_deg * buffer_deg
    for poly in polys:
        if not poly:
            continue
        ring = poly[0]
        for i in range(len(ring) - 1):
            ax, ay = ring[i]
            bx, by = ring[i + 1]
            if _point_to_segment_dist_sq(lon, lat, ax, ay, bx, by) <= max_sq:
                return True
    return False


def is_point_in_au_country(
    country_code: str,
    lat: float,
    lon: float,
    coastal_buffer_deg: float = DEFAULT_COASTAL_BUFFER_DEG,
) -> tuple[bool, str | None]:
    """Check whether WGS-84 (lat, lon) lies inside declared AU country polygon or coastal buffer."""
    code = country_code.strip().upper()
    index = _load_country_polygons()
    if code not in index:
        return False, f"No boundary polygon registered for country code '{code}'."

    polys, (min_lon, min_lat, max_lon, max_lat) = index[code]
    if (
        min_lon - coastal_buffer_deg <= lon <= max_lon + coastal_buffer_deg
        and min_lat - coastal_buffer_deg <= lat <= max_lat + coastal_buffer_deg
    ) and any(_point_in_polygon(lon, lat, poly) for poly in polys):
        return True, None

    for other_code, (other_polys, (omin_x, omin_y, omax_x, omax_y)) in index.items():
        if other_code == code:
            continue
        if (
            omin_x <= lon <= omax_x
            and omin_y <= lat <= omax_y
            and any(_point_in_polygon(lon, lat, poly) for poly in other_polys)
        ):
            return (
                False,
                (
                    f"Coordinates ({lat}, {lon}) fall inside '{other_code}', "
                    f"not declared country '{code}'."
                ),
            )

    if (
        min_lon - coastal_buffer_deg <= lon <= max_lon + coastal_buffer_deg
        and min_lat - coastal_buffer_deg <= lat <= max_lat + coastal_buffer_deg
        and _within_coastal_buffer(lon, lat, polys, coastal_buffer_deg)
    ):
        return True, None

    dist_deg = math.inf
    for poly in polys:
        if not poly:
            continue
        ring = poly[0]
        for i in range(len(ring) - 1):
            ax, ay = ring[i]
            bx, by = ring[i + 1]
            d = math.sqrt(_point_to_segment_dist_sq(lon, lat, ax, ay, bx, by))
            if d < dist_deg:
                dist_deg = d

    return (
        False,
        f"Coordinates ({lat}, {lon}) lie outside '{code}' Natural Earth boundary "
        f"(nearest boundary distance {dist_deg:.2f} deg).",
    )
