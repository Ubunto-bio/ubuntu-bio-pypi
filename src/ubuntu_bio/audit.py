"""External audit checkpoint verification and append-only anchoring client."""

from __future__ import annotations

import base64
import hashlib
import hmac
import json
import time
from pathlib import Path
from typing import Any

import httpx

from ubuntu_bio.config import get_hub_url, get_token
from ubuntu_bio.exceptions import AuthenticationError, UbuntuBioError

# DER DigestInfo prefix for SHA-256 (RFC 8017 Section 9.2)
_SHA256_DIGEST_INFO_PREFIX = bytes.fromhex("3031300d060960864801650304020105000420")


def _b64url_decode(segment: str) -> bytes:
    pad = "=" * (-len(segment) % 4)
    return base64.urlsafe_b64decode(segment + pad)


def _verify_rs256_signature(signing_input: bytes, signature: bytes, jwk: dict[str, Any]) -> bool:
    """Verify RSASSA-PKCS1-v1_5 SHA-256 signature against a JWK public key."""
    if jwk.get("kty") != "RSA":
        return False
    n_b64 = jwk.get("n")
    e_b64 = jwk.get("e")
    if not isinstance(n_b64, str) or not isinstance(e_b64, str):
        return False

    n_bytes = _b64url_decode(n_b64)
    e_bytes = _b64url_decode(e_b64)
    n = int.from_bytes(n_bytes, "big")
    e = int.from_bytes(e_bytes, "big")
    k = (n.bit_length() + 7) // 8
    if k < 64 or len(signature) != k:
        return False

    digest_info = _SHA256_DIGEST_INFO_PREFIX + hashlib.sha256(signing_input).digest()
    ps_len = k - len(digest_info) - 3
    if ps_len < 8:
        return False
    expected_em = b"\x00\x01" + (b"\xff" * ps_len) + b"\x00" + digest_info

    s_int = int.from_bytes(signature, "big")
    m_int = pow(s_int, e, n)
    decoded_em = m_int.to_bytes(k, "big")
    return hmac.compare_digest(decoded_em, expected_em)


def verify_rs256_checkpoint_jwt(
    token: str,
    jwks_data: dict[str, Any],
    *,
    expected_repo: str | None = None,
    now_ts: int | None = None,
) -> dict[str, Any]:
    """Verify an RS256 audit checkpoint JWT against server JWKS keys."""
    parts = token.split(".")
    if len(parts) != 3:
        raise UbuntuBioError("Malformed checkpoint JWT: expected 3 segments.")

    try:
        header = json.loads(_b64url_decode(parts[0]))
        payload: dict[str, Any] = json.loads(_b64url_decode(parts[1]))
        signature = _b64url_decode(parts[2])
    except Exception as exc:
        raise UbuntuBioError(f"Failed decoding checkpoint JWT: {exc}") from exc

    if header.get("alg") != "RS256":
        raise UbuntuBioError(f"Unsupported checkpoint JWT algorithm: {header.get('alg')}")

    keys = jwks_data.get("keys")
    if not isinstance(keys, list) or not keys:
        raise UbuntuBioError("Server JWKS response contained no verification keys.")

    kid = header.get("kid")
    candidates = [k for k in keys if isinstance(k, dict) and (kid is None or k.get("kid") == kid)]
    if not candidates:
        raise UbuntuBioError(f"No JWKS key matched checkpoint kid '{kid}'.")

    signing_input = f"{parts[0]}.{parts[1]}".encode("ascii")
    if not any(_verify_rs256_signature(signing_input, signature, jwk) for jwk in candidates):
        raise UbuntuBioError("Invalid RS256 signature on audit checkpoint JWT.")

    current_ts = int(time.time()) if now_ts is None else now_ts
    exp = payload.get("exp")
    if not isinstance(exp, int | float) or exp < current_ts:
        raise UbuntuBioError("Audit checkpoint JWT has expired.")

    if payload.get("type") != "audit_checkpoint":
        raise UbuntuBioError("Token is not an audit_checkpoint payload.")
    if expected_repo and payload.get("repo") != expected_repo:
        raise UbuntuBioError(
            f"Checkpoint repository scope '{payload.get('repo')}' != '{expected_repo}'."
        )
    if not isinstance(payload.get("entries_checked"), int) or not isinstance(
        payload.get("head_hash"), str
    ):
        raise UbuntuBioError("Audit checkpoint payload missing entries_checked or head_hash.")
    return payload


def load_anchor_ledger(anchor_path: Path) -> list[dict[str, Any]]:
    """Load previously saved external checkpoint anchors from JSONL or JSON."""
    if not anchor_path.exists():
        return []
    text = anchor_path.read_text(encoding="utf-8").strip()
    if not text:
        return []
    if text.startswith("["):
        data = json.loads(text)
        return [item for item in data if isinstance(item, dict)]
    if text.startswith("{") and "\n" not in text:
        single = json.loads(text)
        return [single] if isinstance(single, dict) else []
    records: list[dict[str, Any]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            obj = json.loads(stripped)
            if isinstance(obj, dict):
                records.append(obj)
    return records


def verify_checkpoint_continuity(
    claims: dict[str, Any],
    prior_anchors: list[dict[str, Any]],
    audit_entries: list[dict[str, Any]] | None = None,
) -> list[str]:
    """Detect history truncation or rewriting against independently retained anchors."""
    violations: list[str] = []
    repo_slug = str(claims.get("repo", ""))
    cur_count = int(claims["entries_checked"])
    cur_head = str(claims["head_hash"])

    by_seq: dict[int, dict[str, Any]] = {}
    if audit_entries is not None:
        for row in audit_entries:
            seq = row.get("sequence_no")
            if isinstance(seq, int):
                by_seq[seq] = row

    for anchor in prior_anchors:
        if anchor.get("repo") and anchor["repo"] != repo_slug:
            continue
        prev_count = int(anchor["entries_checked"])
        prev_head = str(anchor["head_hash"])

        if cur_count < prev_count:
            violations.append(
                f"Audit truncation detected: server reported {cur_count} entries "
                f"below retained anchor {prev_count}."
            )
            continue

        if cur_count == prev_count and cur_head != prev_head:
            violations.append(
                f"Audit rewrite detected at sequence {cur_count}: "
                f"server head {cur_head[:12]} != retained anchor {prev_head[:12]}."
            )
            continue

        if audit_entries is not None:
            anchored_row = by_seq.get(prev_count)
            if anchored_row is None:
                violations.append(f"Exported audit log missing anchored sequence_no {prev_count}.")
            elif str(anchored_row.get("entry_hash", "")) != prev_head:
                violations.append(
                    f"Historical audit rewrite detected at anchored sequence {prev_count}: "
                    f"log entry hash {str(anchored_row.get('entry_hash', ''))[:12]} "
                    f"!= retained anchor {prev_head[:12]}."
                )

    return violations


def fetch_and_verify_checkpoint(
    repo_slug: str,
    *,
    save_path: Path | None = None,
    verify_against: Path | None = None,
    token: str | None = None,
    hub_url: str | None = None,
) -> dict[str, Any]:
    """Fetch server checkpoint + JWKS, verify signature and prior anchors, and optionally append."""
    active_token = token or get_token()
    if not active_token:
        raise AuthenticationError("Authentication required. Run 'ubio login' or pass a token.")
    if "/" not in repo_slug:
        raise UbuntuBioError("Repository must be in 'owner/name' format.")

    owner, name = repo_slug.split("/", 1)
    base_url = (hub_url or get_hub_url()).rstrip("/")
    headers = {"Authorization": f"Bearer {active_token}"}

    try:
        jwks_resp = httpx.get(f"{base_url}/.well-known/jwks.json", timeout=15.0)
        cp_resp = httpx.get(
            f"{base_url}/api/repos/{owner}/{name}/audit-log/checkpoint",
            headers=headers,
            timeout=15.0,
        )
    except httpx.RequestError as exc:
        raise UbuntuBioError(f"Failed contacting UbuntuBio Hub: {exc}") from exc

    if jwks_resp.status_code != 200:
        raise UbuntuBioError(f"Failed fetching JWKS ({jwks_resp.status_code}): {jwks_resp.text}")
    if cp_resp.status_code != 200:
        raise UbuntuBioError(f"Checkpoint request failed ({cp_resp.status_code}): {cp_resp.text}")

    cp_data: dict[str, Any] = cp_resp.json()
    claims = verify_rs256_checkpoint_jwt(
        str(cp_data["checkpoint_jwt"]),
        jwks_resp.json(),
        expected_repo=repo_slug,
    )

    anchor_sources: list[dict[str, Any]] = []
    if verify_against is not None:
        anchor_sources.extend(load_anchor_ledger(verify_against))
    if save_path is not None and save_path != verify_against and save_path.exists():
        anchor_sources.extend(load_anchor_ledger(save_path))

    audit_entries: list[dict[str, Any]] | None = None
    if anchor_sources:
        try:
            log_resp = httpx.get(
                f"{base_url}/api/repos/{owner}/{name}/audit-log?format=json",
                headers=headers,
                timeout=30.0,
            )
            if log_resp.status_code == 200 and isinstance(log_resp.json(), list):
                audit_entries = log_resp.json()
        except httpx.RequestError:
            audit_entries = None

    violations = verify_checkpoint_continuity(claims, anchor_sources, audit_entries=audit_entries)
    if violations:
        raise UbuntuBioError("; ".join(violations))

    anchor_record = {
        "repo": repo_slug,
        "entries_checked": claims["entries_checked"],
        "head_hash": claims["head_hash"],
        "kid": cp_data.get("kid"),
        "jti": claims.get("jti"),
        "iat": claims.get("iat"),
        "checkpoint_jwt": cp_data["checkpoint_jwt"],
    }
    if save_path is not None:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with save_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(anchor_record) + "\n")

    return {
        **anchor_record,
        "anchors_checked": len(
            [a for a in anchor_sources if not a.get("repo") or a.get("repo") == repo_slug]
        ),
        "saved_to": str(save_path) if save_path is not None else None,
    }
