"""Biological and schema constants for UbuntuBio client library."""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Final

# 55 African Union (AU) member states (ISO-3166-1 alpha-2, including EH)
AU_COUNTRY_CODES: Final[frozenset[str]] = frozenset(
    {
        "DZ",
        "AO",
        "BJ",
        "BW",
        "BF",
        "BI",
        "CV",
        "CM",
        "CF",
        "TD",
        "KM",
        "CG",
        "CD",
        "CI",
        "DJ",
        "EG",
        "EH",
        "GQ",
        "ER",
        "SZ",
        "ET",
        "GA",
        "GM",
        "GH",
        "GN",
        "GW",
        "KE",
        "LS",
        "LR",
        "LY",
        "MG",
        "MW",
        "ML",
        "MR",
        "MU",
        "MA",
        "MZ",
        "NA",
        "NE",
        "NG",
        "RW",
        "ST",
        "SN",
        "SC",
        "SL",
        "SO",
        "ZA",
        "SS",
        "SD",
        "TZ",
        "TG",
        "TN",
        "UG",
        "ZM",
        "ZW",
    }
)

# Strictly prohibited taxa (human genomic data is forbidden in v1)
PROHIBITED_TAXA: Final[frozenset[int]] = frozenset({9606})

# Permitted biological domains
VALID_DOMAINS: Final[frozenset[str]] = frozenset(
    {
        "bacteria",
        "virus",
        "veterinary",
        "plant",
        "parasite",
        "metagenome",
    }
)

# Permitted access tiers
VALID_ACCESS_TIERS: Final[frozenset[str]] = frozenset(
    {
        "public",
        "gated",
        "private",
    }
)

# Permitted geospatial scope policies and human sequence screening modes
VALID_GEO_POLICIES: Final[frozenset[str]] = frozenset({"au_strict", "allow_outgroups"})
VALID_HUMAN_SCREEN_POLICIES: Final[frozenset[str]] = frozenset({"warn", "strict_reject"})

# Permitted sequencing platforms
VALID_PLATFORMS: Final[frozenset[str]] = frozenset(
    {
        "ILLUMINA",
        "OXFORD_NANOPORE",
        "PACBIO_SMRT",
        "ELEMENT",
        "MGI",
        "ION_TORRENT",
    }
)

# Permitted library strategies
VALID_STRATEGIES: Final[frozenset[str]] = frozenset(
    {
        "WGS",
        "RNA-SEQ",
        "AMPLICON",
        "TARGETED-CAPTURE",
        "METAGENOMIC",
        "OTHER",
    }
)

# Permitted library layouts
VALID_LAYOUTS: Final[frozenset[str]] = frozenset(
    {
        "paired",
        "single",
    }
)

# Permitted library file types
VALID_FILE_TYPES: Final[frozenset[str]] = frozenset(
    {
        "FASTQ",
        "BAM",
        "FASTA",
        "CRAM",
    }
)

# Standard sequence data file extensions for parity checks
SEQUENCE_FILE_EXTENSIONS: Final[tuple[str, ...]] = (
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".fasta",
    ".fasta.gz",
    ".fa",
    ".fa.gz",
    ".bam",
    ".cram",
)

# Mandatory columns for samples.tsv
MANDATORY_SAMPLE_COLUMNS: Final[tuple[str, ...]] = (
    "sample_id",
    "taxon_id",
    "country",
    "collection_date",
    "lat",
    "lon",
    "admin1",
    "host",
    "isolation_source",
)

# Geographic bounding box for the African continent and surrounding islands
AFRICA_LAT_MIN: Final[float] = -35.5
AFRICA_LAT_MAX: Final[float] = 38.0
AFRICA_LON_MIN: Final[float] = -26.0
AFRICA_LON_MAX: Final[float] = 52.0

# Sequence file magic bytes
GZIP_MAGIC: Final[bytes] = b"\x1f\x8b"
BAM_MAGIC: Final[bytes] = b"BAM\x01"
CRAM_MAGIC: Final[bytes] = b"CRAM"

# Human chromosome / assembly header indicators in FASTA/BAM contigs
HUMAN_HEADER_REGEX: Final[re.Pattern[str]] = re.compile(
    r"(?i)\b(homo\s+sapiens|grch3[78]|hg19|hg38|hs37d5|human_g1k)\b|"
    r"^(chr([1-9]|1[0-9]|2[0-2]|[XYM]|MT)|NC_0000(0[1-9]|1[0-9]|2[0-4]))(\.|\b|$)",
    re.IGNORECASE,
)


# Curated strand-specific 25-mer panel tiled across human mtDNA (NC_012920.1),
# 18S rRNA (NR_003286.4), 28S rRNA (NR_003287.4), and AluY (U14574.1) after subtracting
# cross-kingdom bacterial, viral, and protozoan pathogen 25-mers.
def _load_human_kmers() -> frozenset[str]:
    kmer_path = Path(__file__).with_name("human_kmers.txt.gz")
    if kmer_path.exists():
        lines = gzip.decompress(kmer_path.read_bytes()).decode("utf-8").splitlines()
        return frozenset(line.strip().upper() for line in lines if line.strip())
    return frozenset(
        {
            "CGCGACCTCAGATCAGACGTGGCGA",
            "TGCCAGTAGCATATGCTTGTCTCAA",
            "GTAACCCGTTGAACCCCATTCGTGA",
            "AGGAATTGACGGAAGGGCACCACCA",
            "TCCTGCCAGTAGCATATGCTTGTCT",
            "GATCACAGGTCTATCACCCTATTAA",
            "TAGACGGGCTCACATCACCCCATAA",
            "ACCCTATTAACCACTCACGGGAGCT",
        }
    )


HUMAN_SIGNATURE_KMERS: Final[frozenset[str]] = _load_human_kmers()
