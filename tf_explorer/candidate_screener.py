"""
Authoritative Candidate Screening and Purity Exclusion Module for TF-Explorer Phase 4C-3.

Enforces:
1. Pinned cryptographic protection of PHASE4C3_BANNED_EXPOSURE_REGISTER.json.
2. Fail-closed candidate metadata validation (rejecting empty/missing identity,
   missing TF, missing biosample, non-GRCh38 assembly, boolean/fractional coordinates,
   and coordinates exceeding canonical chromosome lengths).
3. Canonical GRCh38 chromosome normalization and coordinate validation.
4. Independent human gene/transcript/TSS annotation verification (via MANE Select & verified transcripts)
   with cryptographic checksum verification of reference data before use.
5. 0-based half-open interval intersection arithmetic (max < min; abutting intervals do NOT overlap).
6. Relational purity exclusion against all 17 authoritative historical records.
7. Genuine ENCODE DCC source authority:
   - Must originate from PRODUCTION_PINNED_SNAPSHOT (cryptographically pinned hash)
     or PRODUCTION_ENCODE_PORTAL (verified HTTPS REST endpoint with response digest).
   - Rejection of self-certified dictionaries, unanchored stores, and synthetic test fixtures fail-closed.
8. Bi-directional matching of biological treatments and genetic modification event sets:
   - Rejection of candidates claiming untreated/wildtype when DCC records perturbation.
   - Rejection of candidates claiming unsupported perturbations when DCC records untreated/wildtype.
9. Required non-empty authoritative file size (> 0) and 32-character hexadecimal MD5 checksum.
10. Independent primary BED byte authentication (when physical BED file is present) vs
    provisional metadata-verified status (when only DCC metadata are available).
11. Collision-safe, atomic ledger identification and writing preventing duplicate SCR-XXXX IDs.
12. Dual cryptographic recording of caller candidate metadata SHA-256 and authoritative DCC response SHA-256.
"""

from __future__ import annotations

import csv
import datetime
import hashlib
import hmac
import json
import os
from pathlib import Path
import re
import secrets
import time
from typing import Any, NamedTuple
import urllib.request
import urllib.error
import urllib.parse
import gzip

_PORTAL_SIGNING_SECRET = secrets.token_hex(32)
_ACQUISITION_TOKEN = secrets.token_hex(32)

def _compute_portal_signature(payload_bytes: bytes) -> str:
    """Computes an internal HMAC signature for verified HTTPS ENCODE Portal responses."""
    return hmac.new(_PORTAL_SIGNING_SECRET.encode("utf-8"), payload_bytes, hashlib.sha256).hexdigest()


FROZEN_BANNED_REGISTER_SHA256 = "7ff48e299bec052a6fe040b133f27fb11408ab10b84e98283e84a270dece20fe"
FROZEN_PRIMARY_PEAK_INDEX_SHA256 = "27228f7041a8b2776668ef92665930fce57122809b6da904d733c12ff6fb4950"
FROZEN_RAW_SNAPSHOTS_SHA256 = "dc7d042e0caf8ac01d4c25788b03d1d50fd0bfe0d33d3ea4415a4c2e53967898"
FROZEN_MANE_SELECT_SHA256 = "f53eb2f55cb4f42b82a806ed8fb27c7140accce0f6c722ac8fc7cf0b62eaef0f"
FROZEN_VERIFIED_TRANSCRIPTS_SHA256 = "395f557ac15e649e6bcf01cb04879502c79c337fa97af697fe93df1466b71c2b"

KNOWN_PINNED_SNAPSHOT_HASHES = {
    FROZEN_PRIMARY_PEAK_INDEX_SHA256,
    FROZEN_RAW_SNAPSHOTS_SHA256,
}

EXPECTED_HISTORICAL_RECORD_IDS = [
    "EXP-4B-01", "EXP-4B-02", "EXP-4B-03", "EXP-4B-04", "EXP-4B-05", "EXP-4B-06",
    "EXP-4C-01", "EXP-4C-02", "EXP-4C-03", "EXP-4C-04", "EXP-4C-05", "EXP-4C-06",
    "EXP-4C-07", "EXP-4C-08", "EXP-4C-09", "EXP-4C-10", "EXP-4C-11"
]

GRCH38_CHROMOSOME_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}

VALID_GRCH38_CHROMOSOMES = set(GRCH38_CHROMOSOME_LENGTHS.keys())

CHROMOSOME_ALIASES = {
    f"{i}": f"chr{i}" for i in range(1, 23)
}
CHROMOSOME_ALIASES.update({
    f"0{i}": f"chr{i}" for i in range(1, 10)
})
CHROMOSOME_ALIASES.update({
    "x": "chrX", "X": "chrX", "chrx": "chrX", "chrX": "chrX", "CHRX": "chrX",
    "y": "chrY", "Y": "chrY", "chry": "chrY", "chrY": "chrY", "CHRY": "chrY",
    "m": "chrM", "M": "chrM", "mt": "chrM", "MT": "chrM", "chrm": "chrM",
    "chrM": "chrM", "chrMT": "chrM", "chrmt": "chrM",
})

APPROVED_PEAK_FORMATS = {
    "narrowpeak", "broadpeak", "bed narrowpeak", "bed broadpeak",
    "optimal idr thresholded peaks", "conservative idr thresholded peaks",
    "pseudoreplicated idr thresholded peaks", "peaks"
}


class AuthorityTrustLevel:
    PRODUCTION_PINNED_SNAPSHOT = "PRODUCTION_PINNED_SNAPSHOT"
    PRODUCTION_ENCODE_PORTAL = "PRODUCTION_ENCODE_PORTAL"
    SYNTHETIC_TEST_FIXTURE = "SYNTHETIC_TEST_FIXTURE"
    UNTRUSTED_CALLER_STORE = "UNTRUSTED_CALLER_STORE"


class RegisterIntegrityError(Exception):
    """Raised when the banned exposure register fails cryptographic or schema verification."""
    pass


class AuthorityIntegrityError(Exception):
    """Raised when authoritative ENCODE source metadata fails integrity or provenance verification."""
    pass


class CandidateValidationError(Exception):
    """Raised when a candidate fails fail-closed input validation."""
    pass


def normalize_chromosome(chrom: Any) -> str | None:
    """
    Normalizes candidate or register chromosome string to canonical GRCh38 format ('chr1'-'chr22', 'chrX', 'chrY', 'chrM').
    Returns None if missing, malformed, or unmapped.
    """
    if chrom is None or not isinstance(chrom, str):
        return None
    raw = chrom.strip()
    if not raw:
        return None

    if raw in VALID_GRCH38_CHROMOSOMES:
        return raw

    if raw in CHROMOSOME_ALIASES:
        return CHROMOSOME_ALIASES[raw]

    lower = raw.lower()
    if lower in CHROMOSOME_ALIASES:
        return CHROMOSOME_ALIASES[lower]

    if lower.startswith("chr"):
        sub = lower[3:]
        if sub in CHROMOSOME_ALIASES:
            return CHROMOSOME_ALIASES[sub]
        if sub.isdigit():
            val = int(sub)
            if 1 <= val <= 22:
                return f"chr{val}"

    return None


def is_valid_integer_coordinate(val: Any) -> bool:
    """
    Validates strictly that a coordinate is an integer coordinate.
    Fails closed on booleans (since in Python bool is an int subclass) and floats/fractions.
    """
    if isinstance(val, bool):
        return False
    if isinstance(val, int):
        return True
    if isinstance(val, str):
        s = val.strip()
        return s.isdigit()
    return False


def intervals_overlap_half_open(start1: int, end1: int, start2: int, end2: int) -> bool:
    """
    Computes whether two 0-based half-open intervals [start1, end1) and [start2, end2) overlap.
    In half-open coordinates:
        max(start1, start2) < min(end1, end2)
    """
    return max(start1, start2) < min(end1, end2)


def compute_candidate_hash(candidate: dict[str, Any]) -> str:
    """Computes deterministic SHA-256 hash of sorted candidate dictionary."""
    serialized = json.dumps(candidate, sort_keys=True, default=str)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def compute_bytes_sha256(data: bytes) -> str:
    """Computes SHA-256 hexadecimal hash of raw bytes."""
    return hashlib.sha256(data).hexdigest()


def compute_bytes_md5(data: bytes) -> str:
    """Computes MD5 hexadecimal hash of raw bytes."""
    return hashlib.md5(data).hexdigest()


def validate_bed_format(bed_bytes: bytes, filename: str = "") -> tuple[bool, str]:
    """
    Validates that file bytes represent a well-formed BED (or gzipped BED) file.
    Verifies gzip header magic bytes and decompression if gzipped,
    and checks that records contain at least 3 tab-separated columns with integer intervals.
    """
    if not bed_bytes or len(bed_bytes) == 0:
        return False, "REJECTED_INVALID_BED_FORMAT: File is empty"

    is_gzip = filename.endswith(".gz") or bed_bytes[:2] == b"\x1f\x8b"
    if is_gzip:
        if bed_bytes[:2] != b"\x1f\x8b":
            return False, "REJECTED_INVALID_BED_FORMAT: File has .gz extension but lacks valid gzip magic header (0x1f 0x8b)"
        try:
            decompressed = gzip.decompress(bed_bytes)
            content = decompressed.decode("utf-8", errors="replace")
        except Exception as e:
            return False, f"REJECTED_INVALID_BED_FORMAT: Gzip decompression failed: {e}"
    else:
        try:
            content = bed_bytes.decode("utf-8", errors="replace")
        except Exception as e:
            return False, f"REJECTED_INVALID_BED_FORMAT: Failed to decode text content: {e}"

    valid_lines = 0
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            return False, f"REJECTED_INVALID_BED_FORMAT: Record has fewer than 3 tab-delimited columns: '{line[:50]}'"
        chrom, start, end = parts[0], parts[1], parts[2]
        if not start.isdigit() or not end.isdigit():
            return False, f"REJECTED_INVALID_BED_FORMAT: Non-integer coordinates in record: [{start}, {end}) on {chrom}"
        if int(start) >= int(end):
            return False, f"REJECTED_INVALID_BED_FORMAT: Inverted or empty interval: [{start}, {end}) on {chrom}"
        valid_lines += 1

    if valid_lines == 0:
        return False, "REJECTED_INVALID_BED_FORMAT: File contains no valid BED genomic interval records"

    return True, "VALID_BED_FORMAT"


# -----------------------------------------------------------------------------
# Gene and Genomic TSS Annotation Reference
# -----------------------------------------------------------------------------

class GeneAnnotationReference:
    """
    Authoritative reference provider for human genes, chromosomes, and canonical TSS coordinates.
    Anchored in MANE Select v1.3 and curated verified transcripts.
    Cryptographically verifies reference file checksums before loading.
    """

    _instance = None

    def __init__(
        self,
        mane_path: Path | str | None = None,
        transcripts_path: Path | str | None = None,
        verify_checksums: bool = True,
    ):
        self.gene_map: dict[str, dict[str, Any]] = {}
        self.verify_checksums = verify_checksums
        self._load_reference_data(mane_path, transcripts_path)

    def _load_reference_data(self, mane_path: Path | str | None, transcripts_path: Path | str | None) -> None:
        root = Path(__file__).resolve().parent.parent.parent

        # 1. Load MANE Select v1.3 (19,288 genes)
        if mane_path is not None:
            m_path = Path(mane_path)
        else:
            m_path = root / "reference_data" / "MANE_v1.3_select.json"

        if not m_path.exists():
            raise FileNotFoundError(f"MANE Select reference file not found at: {m_path}")

        m_bytes = m_path.read_bytes()
        if self.verify_checksums:
            actual_m_sha = hashlib.sha256(m_bytes).hexdigest()
            if actual_m_sha != FROZEN_MANE_SELECT_SHA256:
                raise AuthorityIntegrityError(
                    f"MANE Select reference integrity check failed!\n"
                    f"  Expected: {FROZEN_MANE_SELECT_SHA256}\n"
                    f"  Actual:   {actual_m_sha}"
                )

        data = json.loads(m_bytes.decode("utf-8"))
        for symbol, entry in data.items():
            ref_chr = entry.get("GRCh38_chr", "")
            norm_chr = self._map_ncbi_refseq_chr(ref_chr)
            st = int(entry.get("chr_start", 0))
            en = int(entry.get("chr_end", 0))
            strand = entry.get("chr_strand", "+")
            canonical_tss = st if strand == "+" else en
            self.gene_map[symbol.upper()] = {
                "symbol": symbol.upper(),
                "canonical_chromosome": norm_chr,
                "canonical_tss": canonical_tss,
                "start": st,
                "end": en,
                "strand": strand,
                "hgnc_id": entry.get("HGNC_ID", ""),
                "source": "MANE_v1.3_Select",
            }

        # 2. Supplementary curated verified transcripts (12 genes)
        if transcripts_path is not None:
            t_path = Path(transcripts_path)
        else:
            t_path = root / "PHASE4C1_VERIFIED_TRANSCRIPTS.json"

        if t_path.exists():
            t_bytes = t_path.read_bytes()
            if self.verify_checksums:
                actual_t_sha = hashlib.sha256(t_bytes).hexdigest()
                if actual_t_sha != FROZEN_VERIFIED_TRANSCRIPTS_SHA256:
                    raise AuthorityIntegrityError(
                        f"Verified transcripts reference integrity check failed!\n"
                        f"  Expected: {FROZEN_VERIFIED_TRANSCRIPTS_SHA256}\n"
                        f"  Actual:   {actual_t_sha}"
                    )

            t_data = json.loads(t_bytes.decode("utf-8"))
            for symbol, items in t_data.items():
                if items and isinstance(items, list):
                    first = items[0]
                    sym_upper = symbol.upper()
                    self.gene_map[sym_upper] = {
                        "symbol": sym_upper,
                        "canonical_chromosome": normalize_chromosome(first.get("chromosome")),
                        "canonical_tss": int(first.get("canonical_tss", 0)),
                        "strand": first.get("strand", "+"),
                        "hgnc_id": self.gene_map.get(sym_upper, {}).get("hgnc_id", ""),
                        "source": "PHASE4C1_VERIFIED_TRANSCRIPTS",
                    }

    @staticmethod
    def _map_ncbi_refseq_chr(refseq_acc: str) -> str | None:
        """Maps NCBI RefSeq chromosome accession (e.g. NC_000017.11) to canonical chr (e.g. chr17)."""
        m = re.match(r"^NC_0000(\d\d)\.", refseq_acc)
        if m:
            num = int(m.group(1))
            if 1 <= num <= 22:
                return f"chr{num}"
            if num == 23:
                return "chrX"
            if num == 24:
                return "chrY"
        return None

    def lookup_gene(self, symbol: str) -> dict[str, Any] | None:
        """Looks up canonical gene annotation. Returns None if unknown or unannotated."""
        if not symbol or not isinstance(symbol, str):
            return None
        return self.gene_map.get(symbol.strip().upper())

    def validate_gene_and_coordinates(
        self,
        symbol: str,
        chromosome: str,
        interval_start: int,
        interval_end: int,
        max_promoter_distance_bp: int = 20000,
    ) -> tuple[bool, str, dict[str, Any] | None]:
        """
        Validates gene symbol existence, chromosome match, chromosome length bounds,
        and proximity of evaluated coordinates to the gene's canonical TSS.
        """
        gene_info = self.lookup_gene(symbol)
        if gene_info is None:
            return False, f"REJECTED_UNKNOWN_OR_UNANNOTATED_GENE_SYMBOL: {symbol}", None

        # Chromosome match
        ref_chr = gene_info["canonical_chromosome"]
        if ref_chr != chromosome:
            return False, f"REJECTED_GENE_CHROMOSOME_MISMATCH (candidate: {chromosome}, reference: {ref_chr})", gene_info

        # Chromosome length bounds
        max_chr_len = GRCH38_CHROMOSOME_LENGTHS.get(chromosome, 0)
        if interval_end > max_chr_len:
            return False, f"FAIL_CLOSED_COORDINATES_EXCEED_CHROMOSOME_LENGTH ({interval_end} > {max_chr_len})", gene_info

        # TSS proximity check (promoter must be within promoter window distance of TSS)
        tss = gene_info["canonical_tss"]
        if not (interval_start - max_promoter_distance_bp <= tss <= interval_end + max_promoter_distance_bp):
            return False, (
                f"REJECTED_GENOMIC_COORDINATES_DISCORDANT_WITH_GENE_LOCUS "
                f"(interval [{interval_start}, {interval_end}) is {min(abs(interval_start - tss), abs(interval_end - tss))} bp away from TSS {tss})"
            ), gene_info

        return True, "VALID_GENE_AND_PROMOTER_COORDINATES", gene_info


_GLOBAL_GENE_ANNOTATOR: GeneAnnotationReference | None = None

def get_gene_annotator() -> GeneAnnotationReference:
    global _GLOBAL_GENE_ANNOTATOR
    if _GLOBAL_GENE_ANNOTATOR is None:
        _GLOBAL_GENE_ANNOTATOR = GeneAnnotationReference()
    return _GLOBAL_GENE_ANNOTATOR


# -----------------------------------------------------------------------------
# Encode Provenance Record & Acquisition Pipeline
# -----------------------------------------------------------------------------

class EncodeProvenanceRecord(NamedTuple):
    retrieval_origin: str
    retrieval_timestamp: str
    experiment_accession: str
    file_accession: str
    experiment_status: str
    file_status: str
    dataset: str
    tf_target: str
    biosample_term_name: str
    treatments: list[Any]
    genetic_modifications: list[Any]
    biosample_summary: str
    assembly: str
    file_format: str
    file_format_type: str
    output_type: str
    file_size_bytes: int
    file_md5: str
    raw_response_sha256: str
    provenance_payload: dict[str, Any]
    origin_signature: str = ""
    acquisition_token: str = ""


def fetch_authoritative_encode_provenance(
    experiment_accession: str,
    file_accession: str,
    timeout_seconds: float = 15.0,
) -> EncodeProvenanceRecord:
    """
    Executes an authenticated, audited retrieval from the official ENCODE Portal REST API.
    Enforces HTTPS, official host www.encodeproject.org, HTTP 200, Content-Type application/json,
    computes deterministic response SHA-256 digests, and captures full original metadata
    WITHOUT generating assumed defaults for missing fields.

    Verifies returned experiment/file accessions and mutual dataset membership fail-closed.
    Rejects redirects outside the trusted ENCODE HTTPS origin.
    """
    exp_acc = experiment_accession.strip().upper()
    file_acc = file_accession.strip().upper()

    if not re.match(r"^ENCSR[0-9A-Z]+$", exp_acc):
        raise AuthorityIntegrityError(f"Malformed ENCODE experiment accession: '{exp_acc}'")
    if not re.match(r"^ENCFF[0-9A-Z]+$", file_acc):
        raise AuthorityIntegrityError(f"Malformed ENCODE file accession: '{file_acc}'")

    exp_url = f"https://www.encodeproject.org/experiments/{exp_acc}/?format=json"
    file_url = f"https://www.encodeproject.org/files/{file_acc}/?format=json"

    # Validate request URLs
    for u in (exp_url, file_url):
        p_u = urllib.parse.urlparse(u)
        if p_u.scheme.lower() != "https" or p_u.netloc.lower() not in {"www.encodeproject.org", "encodeproject.org"}:
            raise AuthorityIntegrityError(f"Invalid ENCODE request URL: '{u}'")

    headers = {"Accept": "application/json", "User-Agent": "TF-Explorer-Curator-Authentication-Gate/4.6"}

    def _fetch_endpoint(url: str, acc: str, expected_type: str) -> tuple[bytes, dict[str, Any]]:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=timeout_seconds) as resp:
            # 1. HTTP status code
            status = getattr(resp, "status", getattr(resp, "code", 200))
            if status != 200:
                raise AuthorityIntegrityError(f"ENCODE {expected_type} retrieval failed with HTTP {status} for {url}")

            # 2. Final response URL & redirect validation
            final_url = resp.geturl() if hasattr(resp, "geturl") else url
            p_final = urllib.parse.urlparse(final_url)
            if p_final.scheme.lower() != "https" or p_final.netloc.lower() not in {"www.encodeproject.org", "encodeproject.org"}:
                raise AuthorityIntegrityError(
                    f"Redirect outside trusted ENCODE HTTPS origin rejected: requested '{url}', redirected to '{final_url}'"
                )

            # 3. Content-Type header validation
            content_type = ""
            if hasattr(resp, "headers") and resp.headers:
                if hasattr(resp.headers, "get_content_type"):
                    content_type = resp.headers.get_content_type()
                elif hasattr(resp.headers, "get"):
                    content_type = resp.headers.get("Content-Type", "application/json")
            else:
                content_type = "application/json"

            if "application/json" not in content_type.lower():
                raise AuthorityIntegrityError(
                    f"Invalid Content-Type for ENCODE response: expected 'application/json', got '{content_type}'"
                )

            # 4. Raw response bytes
            raw_b = resp.read()

        try:
            data = json.loads(raw_b.decode("utf-8"))
        except Exception as e:
            raise AuthorityIntegrityError(f"Malformed JSON in ENCODE {expected_type} response from {url}: {e}")

        ret_acc = str(data.get("accession", "")).strip().upper()
        if ret_acc != acc:
            raise AuthorityIntegrityError(
                f"Returned ENCODE {expected_type} accession mismatch: requested '{acc}', received '{ret_acc}'"
            )
        return raw_b, data

    exp_bytes, exp_data = _fetch_endpoint(exp_url, exp_acc, "experiment")
    file_bytes, file_data = _fetch_endpoint(file_url, file_acc, "file")

    # Strict experiment/file dataset relationship verification
    file_dataset = str(file_data.get("dataset", "")).strip()
    expected_dataset = f"/experiments/{exp_acc}/"
    if file_dataset.rstrip("/") != expected_dataset.rstrip("/"):
        raise AuthorityIntegrityError(
            f"Contradictory experiment/file relationship: file {file_acc} belongs to dataset '{file_dataset}', "
            f"expected '{expected_dataset}'"
        )

    # Verify file is referenced in experiment files list if present
    exp_files = exp_data.get("files")
    if isinstance(exp_files, list) and exp_files:
        file_found = any(
            (isinstance(f, dict) and f.get("accession") == file_acc) or
            (isinstance(f, str) and file_acc in f)
            for f in exp_files
        )
        if not file_found:
            raise AuthorityIntegrityError(
                f"Contradictory experiment/file relationship: file {file_acc} not found in experiment {exp_acc} files list"
            )

    combined_raw_bytes = exp_bytes + b"\n---ENCODE_FILE_DELIMITER---\n" + file_bytes
    raw_response_sha256 = hashlib.sha256(combined_raw_bytes).hexdigest()
    origin_signature = _compute_portal_signature(combined_raw_bytes)
    timestamp_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()

    combined_payload = {
        "experiment": exp_data,
        "file": file_data,
        "retrieval_urls": [exp_url, file_url],
        "raw_bytes": combined_raw_bytes,
    }

    # Strict field extraction: do NOT invent defaults for missing values
    exp_status = exp_data.get("status")
    file_status = file_data.get("status")
    if not exp_status or not file_status:
        raise AuthorityIntegrityError("Original ENCODE response lacks mandatory status fields")

    dcc_size = file_data.get("file_size")
    dcc_md5 = file_data.get("md5sum")
    if dcc_size is None or not isinstance(dcc_size, int) or dcc_size <= 0:
        raise AuthorityIntegrityError(f"Original ENCODE file record lacks positive integer file_size: {dcc_size}")
    if not dcc_md5 or not isinstance(dcc_md5, str) or not re.match(r"^[a-fA-F0-9]{32}$", dcc_md5.strip()):
        raise AuthorityIntegrityError(f"Original ENCODE file record lacks valid 32-hex MD5 checksum: {dcc_md5}")

    target_obj = exp_data.get("target") or {}
    tf_label = target_obj.get("label") or target_obj.get("name") or ""
    if not tf_label:
        raise AuthorityIntegrityError(f"Original ENCODE experiment record lacks target factor label")

    biosample_obj = exp_data.get("biosample_ontology") or {}
    bio_term = biosample_obj.get("term_name") or exp_data.get("biosample_summary") or ""
    if not bio_term:
        raise AuthorityIntegrityError(f"Original ENCODE experiment record lacks biosample term name")

    treatments = exp_data.get("treatments") or []
    genetic_mods = exp_data.get("genetic_modifications") or []
    bio_summary = exp_data.get("biosample_summary") or file_data.get("simple_biosample_summary") or ""

    return EncodeProvenanceRecord(
        retrieval_origin=exp_url,
        retrieval_timestamp=timestamp_utc,
        experiment_accession=exp_acc,
        file_accession=file_acc,
        experiment_status=str(exp_status).lower(),
        file_status=str(file_status).lower(),
        dataset=file_dataset,
        tf_target=tf_label.strip(),
        biosample_term_name=bio_term.strip(),
        treatments=treatments,
        genetic_modifications=genetic_mods,
        biosample_summary=bio_summary,
        assembly=str(file_data.get("assembly", "")).strip().upper(),
        file_format=str(file_data.get("file_format", "")).lower(),
        file_format_type=str(file_data.get("file_format_type", "")).lower(),
        output_type=str(file_data.get("output_type", "")).lower(),
        file_size_bytes=dcc_size,
        file_md5=dcc_md5.strip().lower(),
        raw_response_sha256=raw_response_sha256,
        provenance_payload=combined_payload,
        origin_signature=origin_signature,
        acquisition_token=_ACQUISITION_TOKEN,
    )


# -----------------------------------------------------------------------------
# Authoritative ENCODE Metadata Store & Provenance
# -----------------------------------------------------------------------------

class AuthoritativeDCCStore:
    """
    Verifiable source authority provider for ENCODE experiment and file metadata.
    Anchored strictly in:
    1. PRODUCTION_PINNED_SNAPSHOT: cryptographically pinned snapshot files matching frozen hashes.
    2. PRODUCTION_ENCODE_PORTAL: verified HTTPS retrieval records from www.encodeproject.org.
    3. SYNTHETIC_TEST_FIXTURE: explicitly quarantined test fixtures that CANNOT grant production eligibility.

    A supplied bare dictionary or caller-constructed store with arbitrary labels cannot authenticate itself.
    """

    _INGESTION_TOKEN = secrets.token_hex(32)

    def __init__(
        self,
        catalog_records: dict[str, Any],
        provenance_source: str,
        provenance_hash: str,
        trust_level: str = AuthorityTrustLevel.SYNTHETIC_TEST_FIXTURE,
        retrieval_timestamp: str | None = None,
        verified_source_bytes: bytes | None = None,
        internal_signature: str | None = None,
        _internal_token: Any = None,
    ):
        # Strict constructor validation: prevent caller from bypassing source origin checks
        if trust_level == AuthorityTrustLevel.PRODUCTION_ENCODE_PORTAL:
            if _internal_token != self._INGESTION_TOKEN:
                raise AuthorityIntegrityError(
                    "Direct construction with PRODUCTION_ENCODE_PORTAL is prohibited; "
                    "production portal stores must be instantiated via verified HTTPS acquisition"
                )

        elif trust_level == AuthorityTrustLevel.PRODUCTION_PINNED_SNAPSHOT:
            if provenance_hash not in KNOWN_PINNED_SNAPSHOT_HASHES:
                raise AuthorityIntegrityError(
                    f"Cannot assign PRODUCTION_PINNED_SNAPSHOT with unpinned provenance hash: {provenance_hash}\n"
                    f"Must match one of: {KNOWN_PINNED_SNAPSHOT_HASHES}"
                )
            if _internal_token != self._INGESTION_TOKEN:
                raise AuthorityIntegrityError(
                    "Direct construction of PRODUCTION_PINNED_SNAPSHOT with caller-supplied catalog dictionary is prohibited; "
                    "Substituted catalog record rejected: trusted catalog must be derived from authentic pinned source bytes via validated ingestion functions."
                )

        self.catalog = catalog_records
        self.provenance_source = provenance_source
        self.provenance_hash = provenance_hash
        self.trust_level = trust_level
        self.retrieval_timestamp = retrieval_timestamp or datetime.datetime.now(datetime.timezone.utc).isoformat()

    @classmethod
    def from_primary_peak_index(
        cls,
        index_path: Path | str | None = None,
        expected_sha256: str = FROZEN_PRIMARY_PEAK_INDEX_SHA256,
    ) -> AuthoritativeDCCStore:
        """Loads and cryptographically validates reference_data/primary_peak_index.json."""
        root = Path(__file__).resolve().parent.parent.parent
        p = Path(index_path) if index_path else root / "reference_data" / "primary_peak_index.json"
        if not p.exists():
            raise FileNotFoundError(f"Primary peak index not found at: {p}")

        raw_bytes = p.read_bytes()
        actual_sha = hashlib.sha256(raw_bytes).hexdigest()
        if actual_sha != expected_sha256:
            raise AuthorityIntegrityError(
                f"Primary peak index hash mismatch!\nExpected: {expected_sha256}\nActual:   {actual_sha}"
            )

        data = json.loads(raw_bytes.decode("utf-8"))
        catalog: dict[str, Any] = {}
        for file_acc, file_info in data.items():
            if file_acc == "_manifest" or not isinstance(file_info, dict):
                continue
            exp_acc = file_info.get("experiment")
            if not exp_acc:
                continue

            file_meta = {
                "accession": file_acc,
                "dataset": f"/experiments/{exp_acc}/",
                "status": "released",
                "assembly": file_info.get("assembly", "GRCh38"),
                "file_format": "bed",
                "file_format_type": "narrowPeak",
                "file_type": file_info.get("format", "bed narrowPeak"),
                "output_type": "optimal IDR thresholded peaks",
                "file_size": file_info.get("file_size"),
                "md5sum": file_info.get("md5"),
                "sha256": file_info.get("sha256"),
            }

            if exp_acc not in catalog:
                catalog[exp_acc] = {
                    "accession": exp_acc,
                    "status": "released",
                    "target": {"label": file_info.get("target")},
                    "biosample_ontology": {"term_name": file_info.get("biosample")},
                    "files": [file_meta],
                }
            else:
                catalog[exp_acc]["files"].append(file_meta)
            catalog[file_acc] = file_meta

        return cls(
            catalog,
            provenance_source=str(p),
            provenance_hash=actual_sha,
            trust_level=AuthorityTrustLevel.PRODUCTION_PINNED_SNAPSHOT,
            _internal_token=cls._INGESTION_TOKEN,
        )

    @classmethod
    def from_raw_snapshots(
        cls,
        snapshots_path: Path | str | None = None,
        expected_sha256: str = FROZEN_RAW_SNAPSHOTS_SHA256,
    ) -> AuthoritativeDCCStore:
        """Loads and cryptographically validates PHASE4C1_RAW_METADATA_SNAPSHOTS.json."""
        root = Path(__file__).resolve().parent.parent.parent
        p = Path(snapshots_path) if snapshots_path else root / "PHASE4C1_RAW_METADATA_SNAPSHOTS.json"
        if not p.exists():
            raise FileNotFoundError(f"Raw metadata snapshots not found at: {p}")

        raw_bytes = p.read_bytes()
        actual_sha = hashlib.sha256(raw_bytes).hexdigest()
        if actual_sha != expected_sha256:
            raise AuthorityIntegrityError(
                f"Raw metadata snapshots hash mismatch!\nExpected: {expected_sha256}\nActual:   {actual_sha}"
            )

        data = json.loads(raw_bytes.decode("utf-8"))
        catalog: dict[str, Any] = {}
        for file_acc, file_info in data.items():
            if not isinstance(file_info, dict):
                continue
            dataset_ref = file_info.get("dataset", "")
            exp_match = re.search(r"ENCSR[0-9A-Z]+", dataset_ref)
            exp_acc = exp_match.group(0) if exp_match else None
            catalog[file_acc] = file_info
            if exp_acc:
                if exp_acc not in catalog:
                    tf_obj = file_info.get("target")
                    bio_obj = file_info.get("biosample_ontology")
                    catalog[exp_acc] = {
                        "accession": exp_acc,
                        "status": "released",
                        "target": tf_obj,
                        "biosample_ontology": bio_obj,
                        "files": [file_info],
                    }
                else:
                    catalog[exp_acc]["files"].append(file_info)

        return cls(
            catalog,
            provenance_source=str(p),
            provenance_hash=actual_sha,
            trust_level=AuthorityTrustLevel.PRODUCTION_PINNED_SNAPSHOT,
            _internal_token=cls._INGESTION_TOKEN,
        )

    @classmethod
    def from_encode_provenance_record(
        cls,
        record: EncodeProvenanceRecord,
    ) -> AuthoritativeDCCStore:
        """Constructs an authoritative store from a verified ENCODE Portal retrieval record."""
        # 1. Reject caller-created records not produced by verified HTTPS acquisition
        acq_token = getattr(record, "acquisition_token", "")
        if not acq_token or acq_token != _ACQUISITION_TOKEN:
            raise AuthorityIntegrityError(
                "Caller-created EncodeProvenanceRecord cannot grant production authority: "
                "record was not created by validated HTTPS acquisition pipeline"
            )

        raw_bytes = record.provenance_payload.get("raw_bytes") if isinstance(record.provenance_payload, dict) else None
        if not raw_bytes or not isinstance(raw_bytes, (bytes, bytearray)):
            raise AuthorityIntegrityError(
                "EncodeProvenanceRecord lacks verified raw response bytes; cannot grant production authority"
            )

        if b"\n---ENCODE_FILE_DELIMITER---\n" not in raw_bytes:
            raise AuthorityIntegrityError("Malformed raw response bytes in EncodeProvenanceRecord")

        exp_bytes, file_bytes = raw_bytes.split(b"\n---ENCODE_FILE_DELIMITER---\n", 1)
        try:
            exp_data = json.loads(exp_bytes.decode("utf-8"))
            file_data = json.loads(file_bytes.decode("utf-8"))
        except Exception as e:
            raise AuthorityIntegrityError(f"Failed to parse original ENCODE response bytes: {e}")

        # 2. Check for mutation of record identity fields against parsed source bytes
        parsed_exp_acc = str(exp_data.get("accession", "")).strip().upper()
        parsed_file_acc = str(file_data.get("accession", "")).strip().upper()
        if parsed_exp_acc != record.experiment_accession.strip().upper():
            raise AuthorityIntegrityError(
                f"Mutation of signed-record identity fields detected: experiment accession '{record.experiment_accession}' "
                f"contradicts original source bytes '{parsed_exp_acc}'"
            )
        if parsed_file_acc != record.file_accession.strip().upper():
            raise AuthorityIntegrityError(
                f"Mutation of signed-record identity fields detected: file accession '{record.file_accession}' "
                f"contradicts original source bytes '{parsed_file_acc}'"
            )

        target_obj = exp_data.get("target") or {}
        parsed_tf = target_obj.get("label") or target_obj.get("name") or ""
        if parsed_tf.strip().upper() != record.tf_target.strip().upper():
            raise AuthorityIntegrityError(
                f"Mutation of signed-record identity fields detected: record tf_target '{record.tf_target}' "
                f"contradicts original source bytes '{parsed_tf}'"
            )

        parsed_size = file_data.get("file_size")
        if parsed_size != record.file_size_bytes:
            raise AuthorityIntegrityError(
                f"Mutation of signed-record identity fields detected: record file_size_bytes '{record.file_size_bytes}' "
                f"contradicts original source bytes '{parsed_size}'"
            )

        parsed_md5 = str(file_data.get("md5sum", "")).strip().lower()
        if parsed_md5 != record.file_md5.strip().lower():
            raise AuthorityIntegrityError(
                f"Mutation of signed-record identity fields detected: record file_md5 '{record.file_md5}' "
                f"contradicts original source bytes '{parsed_md5}'"
            )

        # 3. Enforce experiment/file relationship
        expected_ds = f"/experiments/{parsed_exp_acc}/"
        rec_ds = str(file_data.get("dataset", "")).strip()
        if rec_ds.rstrip("/") != expected_ds.rstrip("/"):
            raise AuthorityIntegrityError(
                f"Contradictory experiment/file relationship: file {parsed_file_acc} belongs to dataset '{rec_ds}', "
                f"expected '{expected_ds}'"
            )

        # 4. Build catalog bound directly to parsed source records
        file_meta = {
            "accession": parsed_file_acc,
            "dataset": rec_ds,
            "status": str(file_data.get("status", "")).lower(),
            "assembly": str(file_data.get("assembly", "")).strip().upper(),
            "file_format": str(file_data.get("file_format", "")).lower(),
            "file_format_type": str(file_data.get("file_format_type", "")).lower(),
            "file_type": str(file_data.get("file_type", file_data.get("file_format", ""))).lower(),
            "output_type": str(file_data.get("output_type", "")).lower(),
            "file_size": parsed_size,
            "md5sum": parsed_md5,
        }
        exp_meta = {
            "accession": parsed_exp_acc,
            "status": str(exp_data.get("status", "")).lower(),
            "target": {"label": parsed_tf},
            "biosample_ontology": exp_data.get("biosample_ontology") or {},
            "biosample_summary": exp_data.get("biosample_summary") or "",
            "treatments": exp_data.get("treatments") or [],
            "genetic_modifications": exp_data.get("genetic_modifications") or [],
            "files": [file_meta],
        }
        catalog = {
            parsed_exp_acc: exp_meta,
            parsed_file_acc: file_meta,
        }
        calc_hash = hashlib.sha256(raw_bytes).hexdigest()
        return cls(
            catalog_records=catalog,
            provenance_source=record.retrieval_origin,
            provenance_hash=calc_hash,
            trust_level=AuthorityTrustLevel.PRODUCTION_ENCODE_PORTAL,
            retrieval_timestamp=record.retrieval_timestamp,
            _internal_token=cls._INGESTION_TOKEN,
        )

    @classmethod
    def from_verified_fixture(
        cls,
        fixture_records: dict[str, Any],
        fixture_identifier: str = "INDEPENDENT_TEST_FIXTURE",
    ) -> AuthoritativeDCCStore:
        """
        Constructs an authority store from synthetic test fixture records.
        Strictly labeled SYNTHETIC_TEST_FIXTURE; cannot grant production eligibility.
        """
        serialized = json.dumps(fixture_records, sort_keys=True, default=str).encode("utf-8")
        h = hashlib.sha256(serialized).hexdigest()
        return cls(
            catalog_records=fixture_records,
            provenance_source=f"SYNTHETIC_TEST_FIXTURE:{fixture_identifier}",
            provenance_hash=h,
            trust_level=AuthorityTrustLevel.SYNTHETIC_TEST_FIXTURE,
        )

    def get_experiment_and_file(
        self,
        exp_accession: str,
        file_accession: str,
    ) -> tuple[dict[str, Any] | None, dict[str, Any] | None, str]:
        """
        Retrieves experiment and file records along with deterministic response hash.
        Returns (exp_record, file_record, authority_response_hash).
        """
        exp_acc = exp_accession.strip().upper()
        file_acc = file_accession.strip().upper()

        exp_record = self.catalog.get(exp_acc)
        if exp_record is None:
            return None, None, "NONE_EXPERIMENT_NOT_FOUND"

        file_record = None
        files_list = exp_record.get("files", [])
        if isinstance(files_list, list):
            for f in files_list:
                if isinstance(f, dict) and f.get("accession") == file_acc:
                    file_record = f
                    break

        if file_record is None and file_acc in self.catalog:
            direct_f = self.catalog[file_acc]
            if isinstance(direct_f, dict):
                dataset = str(direct_f.get("dataset", ""))
                if exp_acc in dataset:
                    file_record = direct_f

        if file_record is None:
            return exp_record, None, "NONE_FILE_NOT_ASSOCIATED"

        resp_payload = {
            "experiment": exp_record,
            "file": file_record,
            "provenance_source": self.provenance_source,
            "provenance_hash": self.provenance_hash,
            "trust_level": self.trust_level,
        }
        resp_hash = hashlib.sha256(json.dumps(resp_payload, sort_keys=True, default=str).encode("utf-8")).hexdigest()

        return exp_record, file_record, resp_hash


# -----------------------------------------------------------------------------
# Cross-Platform Collision-Safe File Lock for Curator Ledger
# -----------------------------------------------------------------------------

class LedgerFileLock:
    """
    Cross-platform atomic file lock using os.open(..., O_CREAT | O_EXCL).
    Guarantees that concurrent curator instances cannot collide on sequential screening IDs.
    """

    def __init__(self, ledger_path: Path, timeout_seconds: float = 10.0):
        self.lock_path = ledger_path.parent / f"{ledger_path.name}.lock"
        self.timeout = timeout_seconds
        self.fd: int | None = None

    def __enter__(self):
        start = time.time()
        while True:
            try:
                self.fd = os.open(str(self.lock_path), os.O_CREAT | os.O_EXCL | os.O_RDWR)
                return self
            except FileExistsError:
                if time.time() - start > self.timeout:
                    try:
                        mtime = self.lock_path.stat().st_mtime
                        if time.time() - mtime > 60.0:
                            self.lock_path.unlink(missing_ok=True)
                            continue
                    except Exception:
                        pass
                    raise TimeoutError(f"Could not acquire ledger lock on {self.lock_path} within {self.timeout}s")
                time.sleep(0.005)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.fd is not None:
            try:
                os.close(self.fd)
            except Exception:
                pass
            try:
                self.lock_path.unlink(missing_ok=True)
            except Exception:
                pass


# -----------------------------------------------------------------------------
# Candidate Metadata Validation (Fail-Closed)
# -----------------------------------------------------------------------------

def validate_candidate_metadata(candidate: Any) -> tuple[bool, str, dict[str, Any] | None]:
    """
    Validates that a candidate record has non-empty, well-formed identity metadata.
    Fails closed on missing or malformed candidate identity.
    Returns (is_valid, validation_code, normalized_dict).
    """
    if not isinstance(candidate, dict) or not candidate:
        return False, "FAIL_CLOSED_EMPTY_OR_NON_DICT_CANDIDATE", None

    # Required string fields
    exp_acc = candidate.get("experiment_accession")
    if not isinstance(exp_acc, str) or not exp_acc.strip():
        return False, "FAIL_CLOSED_MISSING_OR_EMPTY_EXPERIMENT_ACCESSION", None

    file_acc = candidate.get("file_accession")
    if not isinstance(file_acc, str) or not file_acc.strip():
        return False, "FAIL_CLOSED_MISSING_OR_EMPTY_FILE_ACCESSION", None

    gene = candidate.get("gene_symbol")
    if not isinstance(gene, str) or not gene.strip():
        return False, "FAIL_CLOSED_MISSING_OR_EMPTY_GENE_SYMBOL", None

    # Transcription Factor (mandatory)
    tf_target = candidate.get("tf_target") or candidate.get("transcription_factor")
    if not isinstance(tf_target, str) or not tf_target.strip():
        return False, "FAIL_CLOSED_MISSING_OR_EMPTY_TF_TARGET", None

    # Biosample (mandatory)
    biosample = candidate.get("biosample") or candidate.get("biosample_lineage")
    if not isinstance(biosample, str) or not biosample.strip():
        return False, "FAIL_CLOSED_MISSING_OR_EMPTY_BIOSAMPLE", None

    # Assembly (mandatory, must be GRCh38)
    assembly = candidate.get("assembly")
    if not isinstance(assembly, str) or assembly.strip().upper() != "GRCH38":
        return False, "FAIL_CLOSED_INVALID_OR_MISSING_ASSEMBLY", None

    # Chromosome normalization
    chrom_raw = candidate.get("chromosome")
    norm_chrom = normalize_chromosome(chrom_raw)
    if norm_chrom is None:
        return False, "FAIL_CLOSED_INVALID_OR_MISSING_CHROMOSOME", None

    # Coordinate interval validation
    start_raw = candidate.get("interval_start")
    end_raw = candidate.get("interval_end")

    if start_raw is None or end_raw is None:
        return False, "FAIL_CLOSED_MISSING_INTERVAL_COORDINATES", None

    if not is_valid_integer_coordinate(start_raw) or not is_valid_integer_coordinate(end_raw):
        return False, "FAIL_CLOSED_NON_INTEGER_OR_BOOLEAN_COORDINATES", None

    cand_start = int(start_raw)
    cand_end = int(end_raw)

    if cand_start < 0 or cand_end < 0:
        return False, "FAIL_CLOSED_NEGATIVE_COORDINATES", None

    if cand_start >= cand_end:
        return False, "FAIL_CLOSED_INVERTED_OR_EMPTY_INTERVAL", None

    max_chr_len = GRCH38_CHROMOSOME_LENGTHS.get(norm_chrom, 0)
    if cand_end > max_chr_len:
        return False, "FAIL_CLOSED_COORDINATES_EXCEED_CHROMOSOME_LENGTH", None

    norm_cand = dict(candidate)
    norm_cand["experiment_accession"] = exp_acc.strip().upper()
    norm_cand["file_accession"] = file_acc.strip().upper()
    norm_cand["gene_symbol"] = gene.strip().upper()
    norm_cand["tf_target"] = tf_target.strip().upper()
    norm_cand["biosample"] = biosample.strip()
    norm_cand["assembly"] = "GRCh38"
    norm_cand["chromosome"] = norm_chrom
    norm_cand["interval_start"] = cand_start
    norm_cand["interval_end"] = cand_end

    return True, "VALID_CANDIDATE_METADATA", norm_cand


# -----------------------------------------------------------------------------
# Screening Result Representation
# -----------------------------------------------------------------------------

class ScreeningResult(NamedTuple):
    is_admissible: bool
    evaluation_decision: str  # AUTHENTICATED_ELIGIBLE_CANDIDATE, PROVISIONAL_METADATA_VERIFIED_CANDIDATE, UNEXPOSED_BY_REGISTER_CHECK, REJECTED_*
    reason_code: str
    matched_record_id: str | None
    sha256_candidate_metadata: str
    sha256_authority_response: str
    screening_id: str
    candidate_summary: dict[str, Any]


# -----------------------------------------------------------------------------
# Banned Exposure Register Loading & Leakage Evaluation
# -----------------------------------------------------------------------------

def load_authoritative_register(
    register_path: str | Path | None = None,
    verify_frozen_hash: bool = True,
) -> list[dict[str, Any]]:
    """
    Loads and cryptographically validates the authoritative banned exposure register.
    Fails closed if the file is missing, altered, or tampered with.
    """
    if register_path is not None:
        p = Path(register_path)
        if not p.exists():
            raise FileNotFoundError(f"Explicitly specified banned register path does not exist: {register_path}")
    else:
        p = Path("PHASE4C3_BANNED_EXPOSURE_REGISTER.json")
        if not p.exists():
            p = Path(__file__).resolve().parent.parent.parent / "PHASE4C3_BANNED_EXPOSURE_REGISTER.json"
        if not p.exists():
            raise FileNotFoundError(f"Default banned exposure register not found at: {p}")

    with open(p, "rb") as f:
        raw_bytes = f.read()

    computed_sha = hashlib.sha256(raw_bytes).hexdigest()
    if verify_frozen_hash and computed_sha != FROZEN_BANNED_REGISTER_SHA256:
        raise RegisterIntegrityError(
            f"Banned register integrity verification failed!\n"
            f"  Computed SHA-256: {computed_sha}\n"
            f"  Frozen SHA-256:   {FROZEN_BANNED_REGISTER_SHA256}"
        )

    try:
        data = json.loads(raw_bytes.decode("utf-8"))
    except Exception as e:
        raise RegisterIntegrityError(f"Malformed register JSON: {e}") from e

    records = data.get("banned_records", [])
    if len(records) != 17:
        raise RegisterIntegrityError(
            f"Banned exposure register must contain exactly 17 authoritative historical records, found {len(records)}"
        )

    loaded_ids = [r.get("record_id") for r in records]
    if loaded_ids != EXPECTED_HISTORICAL_RECORD_IDS:
        raise RegisterIntegrityError(
            f"Banned exposure register record IDs do not match expected historical sequence!\n"
            f"  Loaded:   {loaded_ids}\n"
            f"  Expected: {EXPECTED_HISTORICAL_RECORD_IDS}"
        )

    for r in records:
        for req_field in ["record_id", "experiment_accession", "file_accession", "gene_symbol", "chromosome", "evaluated_promoter_interval"]:
            if req_field not in r:
                raise RegisterIntegrityError(f"Record {r.get('record_id')} missing required field: {req_field}")
        interval = r["evaluated_promoter_interval"]
        if "start" not in interval or "end" not in interval:
            raise RegisterIntegrityError(f"Record {r.get('record_id')} missing start/end in evaluated_promoter_interval")

    return records


def evaluate_candidate_for_leakage(
    candidate: dict[str, Any],
    register_records: list[dict[str, Any]] | None = None,
    verify_frozen_hash: bool = True,
) -> tuple[bool, str, str | None]:
    """
    Evaluates candidate record against the authoritative exposure register.
    Enforces fail-closed metadata validation prior to unexposed acceptance.
    Returns (is_unexposed, reason_code, matched_record_id).
    """
    is_valid, val_code, norm_cand = validate_candidate_metadata(candidate)
    if not is_valid or norm_cand is None:
        return False, val_code, None

    if register_records is None:
        register_records = load_authoritative_register(verify_frozen_hash=verify_frozen_hash)

    cand_exp = norm_cand["experiment_accession"]
    cand_file = norm_cand["file_accession"]
    cand_gene = norm_cand["gene_symbol"]
    cand_hgnc = str(norm_cand.get("hgnc_id", "")).strip().upper()
    cand_chr = norm_cand["chromosome"]
    cand_start = norm_cand["interval_start"]
    cand_end = norm_cand["interval_end"]

    for rec in register_records:
        rec_id = rec["record_id"]

        # 1. Direct Experiment Accession match
        if cand_exp == rec["experiment_accession"].upper():
            return False, "EXCLUDED_BANNED_EXPERIMENT_ACCESSION", rec_id

        # 2. Direct File Accession match
        if cand_file == rec["file_accession"].upper():
            return False, "EXCLUDED_BANNED_FILE_ACCESSION", rec_id

        # 3. Direct Gene Symbol match
        if cand_gene == rec["gene_symbol"].upper():
            return False, "EXCLUDED_BANNED_GENE_IDENTITY", rec_id

        # 4. Direct HGNC ID match
        rec_hgnc = str(rec.get("hgnc_id", "")).strip().upper()
        if cand_hgnc and rec_hgnc and cand_hgnc == rec_hgnc:
            return False, "EXCLUDED_BANNED_HGNC_ID", rec_id

        # 5. Genomic Coordinate Interval Overlap (half-open: max < min)
        rec_chr = normalize_chromosome(rec.get("chromosome"))
        rec_interval = rec.get("evaluated_promoter_interval", {})
        rec_start = rec_interval.get("start")
        rec_end = rec_interval.get("end")

        if rec_chr and rec_chr == cand_chr and rec_start is not None and rec_end is not None:
            if intervals_overlap_half_open(cand_start, cand_end, rec_start, rec_end):
                return False, "EXCLUDED_COORDINATE_INTERVAL_OVERLAP", rec_id

    return True, "UNEXPOSED_BY_REGISTER_CHECK", None


# -----------------------------------------------------------------------------
# Biological Perturbation Normalization & Comparison Helpers
# -----------------------------------------------------------------------------

def normalize_treatment_term(term: Any) -> str:
    """Normalizes treatment term by lowercasing, removing concentration/duration, and stripping punctuation."""
    if not term or not isinstance(term, str):
        return ""
    t = term.lower().strip()
    t = re.sub(r"\b\d+(\.\d+)?\s*(nm|um|mm|m|ug/ml|mg/ml|ng/ml|pm)\b", "", t)
    t = re.sub(r"\bfor\s+\d+(\.\d+)?\s*(hour|hr|h|min|minute|day|d)s?\b", "", t)
    t = re.sub(r"[,\-_/]", " ", t)
    return " ".join(t.split())


def extract_treatment_terms(exp_record: dict[str, Any]) -> set[str]:
    """Extracts all normalized treatment terms from DCC experiment record."""
    terms: set[str] = set()
    treatments = exp_record.get("treatments") or []
    if isinstance(treatments, list):
        for item in treatments:
            if isinstance(item, dict):
                for k in ("treatment_term_name", "treatment_type", "description"):
                    v = item.get(k)
                    if v and isinstance(v, str):
                        norm = normalize_treatment_term(v)
                        if norm:
                            terms.add(norm)
            elif isinstance(item, str):
                norm = normalize_treatment_term(item)
                if norm:
                    terms.add(norm)
    summary = (str(exp_record.get("biosample_summary", "")) + " " + str(exp_record.get("simple_biosample_summary", ""))).lower()
    for drug in ("dexamethasone", "doxorubicin", "tamoxifen", "estradiol", "etoposide", "dmso", "interferon", "tgf-beta", "ethanol"):
        if drug in summary:
            terms.add(drug)
    return terms


def extract_candidate_treatment_terms(norm_cand: dict[str, Any]) -> set[str]:
    """Extracts normalized treatment terms from candidate metadata."""
    terms: set[str] = set()
    raw = norm_cand.get("treatment_context") or norm_cand.get("treatments") or norm_cand.get("treatment")
    if isinstance(raw, str):
        cleaned = raw.strip()
        lower_c = cleaned.lower()
        if any(kw in lower_c for kw in ("none", "untreated", "control", "vehicle", "wild-type", "unmodified")) or not cleaned:
            return set()
        norm = normalize_treatment_term(cleaned)
        if norm:
            terms.add(norm)
    elif isinstance(raw, list):
        for item in raw:
            if isinstance(item, str):
                cleaned = item.strip()
                lower_c = cleaned.lower()
                if any(kw in lower_c for kw in ("none", "untreated", "control", "vehicle", "wild-type", "unmodified")) or not cleaned:
                    continue
                norm = normalize_treatment_term(cleaned)
                if norm:
                    terms.add(norm)
            elif isinstance(item, dict):
                for k in ("treatment_term_name", "treatment_type", "description"):
                    v = item.get(k)
                    if v and isinstance(v, str):
                        cleaned = v.strip()
                        lower_c = cleaned.lower()
                        if any(kw in lower_c for kw in ("none", "untreated", "control", "vehicle", "wild-type", "unmodified")) or not cleaned:
                            continue
                        norm = normalize_treatment_term(cleaned)
                        if norm:
                            terms.add(norm)
    return terms


def _extract_mod_components(text: str) -> tuple[str, str]:
    """Extracts (method, target_factor) from genetic modification string."""
    lower = text.lower()
    method = ""
    for m in ("crispr", "shrna", "sirna", "talen", "zinc-finger", "tagging"):
        if m in lower:
            method = m
            break
    words = re.findall(r"\b[a-zA-Z0-9]{2,10}\b", lower)
    skip = {
        "crispr", "shrna", "sirna", "talen", "knockout", "knockdown", "deletion", "insertion",
        "tagging", "with", "and", "the", "for", "cell", "treated", "untreated", "wild", "type",
        "human", "mcf7", "k562", "hepg2", "a549", "hct116", "gm12878", "none"
    }
    target = ""
    for w in words:
        if w not in skip and not w.isdigit():
            target = w
            break
    return method, target


def extract_genetic_modification_tuples(exp_record: dict[str, Any]) -> set[tuple[str, str]]:
    """Extracts normalized (method, target) tuples from DCC experiment record."""
    mods: set[tuple[str, str]] = set()
    gm_list = exp_record.get("genetic_modifications") or []
    if isinstance(gm_list, list):
        for item in gm_list:
            if isinstance(item, dict):
                method = str(item.get("method", "")).strip().lower()
                target = str(item.get("target") or item.get("modified_site_by_target_id") or "").strip().lower()
                if not method or not target:
                    sub_m, sub_t = _extract_mod_components(str(item))
                    method = method or sub_m
                    target = target or sub_t
                if method or target:
                    mods.add((method, target))
            elif isinstance(item, str):
                m, t = _extract_mod_components(item)
                if m or t:
                    mods.add((m, t))
    summary = (str(exp_record.get("biosample_summary", "")) + " " + str(exp_record.get("simple_biosample_summary", ""))).lower()
    if not mods and any(kw in summary for kw in ("crispr", "shrna", "knockout", "knockdown", "deletion")):
        m, t = _extract_mod_components(summary)
        if m or t:
            mods.add((m, t))
    return mods


def extract_candidate_genetic_modification_tuples(norm_cand: dict[str, Any]) -> set[tuple[str, str]]:
    """Extracts normalized (method, target) tuples from candidate metadata."""
    mods: set[tuple[str, str]] = set()
    raw = norm_cand.get("genetic_modification_context") or norm_cand.get("genetic_modifications") or norm_cand.get("genetic_modification")
    if isinstance(raw, str):
        cleaned = raw.strip()
        lower_c = cleaned.lower()
        if any(kw in lower_c for kw in ("none", "wild-type", "unmodified", "control", "untreated")) or not cleaned:
            return set()
        m, t = _extract_mod_components(cleaned)
        if m or t:
            mods.add((m, t))
    elif isinstance(raw, list):
        for item in raw:
            if isinstance(item, str):
                cleaned = item.strip()
                lower_c = cleaned.lower()
                if any(kw in lower_c for kw in ("none", "wild-type", "unmodified", "control", "untreated")) or not cleaned:
                    continue
                m, t = _extract_mod_components(cleaned)
                if m or t:
                    mods.add((m, t))
            elif isinstance(item, dict):
                method = str(item.get("method", "")).strip().lower()
                target = str(item.get("target") or item.get("modified_site_by_target_id") or "").strip().lower()
                if not method or not target:
                    sub_m, sub_t = _extract_mod_components(str(item))
                    method = method or sub_m
                    target = target or sub_t
                if method or target:
                    mods.add((method, target))
    return mods


# -----------------------------------------------------------------------------
# Authoritative Source Verification (Genuine ENCODE Provenance)
# -----------------------------------------------------------------------------

def verify_authoritative_encode_identity(
    norm_cand: dict[str, Any],
    authority_store: AuthoritativeDCCStore | Any,
    gene_annotator: GeneAnnotationReference | None = None,
) -> tuple[bool, str, dict[str, Any] | None, str]:
    """
    Verifies candidate metadata against authoritative DCC sources with verifiable provenance.
    Enforces:
    1. Rejection of unanchored/bare dictionaries fail-closed.
    2. Experiment & file membership, dataset relationship, and released status.
    3. Approved peak file format and GRCh38 assembly.
    4. Required non-empty authoritative file size (> 0) and 32-hex MD5.
    5. Exact TF target concordance.
    6. Exact biosample concordance (no substring matching).
    7. Exact biological perturbation identity matching (treatments and genetic modifications).
    8. Bi-directional matching of biological treatments and genetic modifications:
       - Rejects untreated/wildtype claims when DCC records perturbations.
       - Rejects unsupported perturbation claims when DCC records untreated/wildtype baseline.
    9. Independent human gene/transcript/TSS annotation verification.

    Returns:
        (is_valid, validation_code, verified_details, authority_response_hash)
    """
    if not isinstance(authority_store, AuthoritativeDCCStore):
        return False, "REJECTED_UNTRUSTED_AUTHORITY_SOURCE: Supplied dictionary has no verifiable provenance anchor", None, "NONE_UNAUTHENTICATED"

    exp_acc = norm_cand["experiment_accession"]
    file_acc = norm_cand["file_accession"]

    # Retrieve authoritative experiment and file records
    exp_record, file_record, auth_hash = authority_store.get_experiment_and_file(exp_acc, file_acc)
    if exp_record is None:
        return False, f"REJECTED_EXPERIMENT_NOT_FOUND_IN_AUTHORITY: {exp_acc}", None, auth_hash
    if file_record is None:
        return False, f"REJECTED_FILE_NOT_ASSOCIATED_WITH_EXPERIMENT ({file_acc} not in {exp_acc})", None, auth_hash

    # Experiment/file dataset relationship verification
    file_dataset = str(file_record.get("dataset", "")).strip()
    if exp_acc not in file_dataset:
        return False, (
            f"REJECTED_FILE_NOT_ASSOCIATED_WITH_EXPERIMENT "
            f"(file {file_acc} belongs to dataset '{file_dataset}', expected '/experiments/{exp_acc}/')"
        ), None, auth_hash

    # Experiment released status
    exp_status = str(exp_record.get("status", "")).lower()
    if exp_status != "released":
        return False, f"REJECTED_REVOKED_OR_NON_RELEASED_EXPERIMENT_STATUS: {exp_status}", None, auth_hash

    # File released status
    file_status = str(file_record.get("status", "")).lower()
    if file_status != "released":
        return False, f"REJECTED_REVOKED_OR_NON_RELEASED_FILE_STATUS: {file_status}", None, auth_hash

    # File format and peak type verification
    file_fmt = str(file_record.get("file_format", "")).lower()
    file_type = str(file_record.get("file_type", "")).lower()
    output_type = str(file_record.get("output_type", "")).lower()
    fmt_type = str(file_record.get("file_format_type", "")).lower()

    if file_fmt != "bed":
        return False, f"REJECTED_NON_BED_FILE_FORMAT: {file_fmt} (required: bed)", None, auth_hash

    is_approved_peak = any(
        ap in file_type or ap in output_type or ap in fmt_type
        for ap in APPROVED_PEAK_FORMATS
    )
    if not is_approved_peak:
        return False, f"REJECTED_UNAPPROVED_FILE_FORMAT: file_type='{file_type}', output_type='{output_type}' (required approved peak format)", None, auth_hash

    # Assembly verification
    file_assembly = str(file_record.get("assembly", "")).strip().upper()
    if file_assembly != "GRCH38":
        return False, f"REJECTED_INCORRECT_ASSEMBLY (file assembly: {file_assembly}, required: GRCh38)", None, auth_hash

    # Non-empty authoritative file size and MD5 checksum
    dcc_size = file_record.get("file_size")
    dcc_md5 = file_record.get("md5sum")

    if dcc_size is None or not isinstance(dcc_size, int) or dcc_size <= 0:
        return False, f"REJECTED_MISSING_OR_INVALID_AUTHORITATIVE_FILE_SIZE (size: {dcc_size})", None, auth_hash

    if not dcc_md5 or not isinstance(dcc_md5, str) or not re.match(r"^[a-fA-F0-9]{32}$", dcc_md5.strip()):
        return False, f"REJECTED_MISSING_OR_INVALID_AUTHORITATIVE_MD5 (md5: {dcc_md5})", None, auth_hash

    # Candidate file size match (if specified)
    cand_size = norm_cand.get("file_size_bytes") or norm_cand.get("file_size")
    if cand_size is not None:
        if not is_valid_integer_coordinate(cand_size) or int(cand_size) != dcc_size:
            return False, f"REJECTED_FILE_SIZE_MISMATCH (candidate: {cand_size}, DCC: {dcc_size})", None, auth_hash

    # Candidate MD5 match (if specified)
    cand_md5 = norm_cand.get("authoritative_md5") or norm_cand.get("md5sum")
    if cand_md5 is not None:
        if str(cand_md5).strip().lower() != dcc_md5.strip().lower():
            return False, f"REJECTED_CHECKSUM_MISMATCH (candidate: {cand_md5}, DCC: {dcc_md5})", None, auth_hash

    # Exact TF target concordance
    dcc_tf = ""
    target_obj = exp_record.get("target")
    if isinstance(target_obj, dict):
        dcc_tf = target_obj.get("label") or target_obj.get("name") or ""
    elif isinstance(target_obj, str):
        dcc_tf = target_obj

    if not dcc_tf or dcc_tf.strip().upper() != norm_cand["tf_target"].upper():
        return False, f"REJECTED_CONTRADICTORY_TF_TARGET (candidate: {norm_cand['tf_target']}, DCC authority: {dcc_tf})", None, auth_hash

    # Exact biosample concordance (strictly no substring matching)
    dcc_biosample = ""
    biosample_ont = exp_record.get("biosample_ontology")
    if isinstance(biosample_ont, dict):
        dcc_biosample = biosample_ont.get("term_name") or ""
    if not dcc_biosample:
        dcc_biosample = exp_record.get("biosample_summary") or exp_record.get("biosample") or ""

    cand_biosample_clean = norm_cand["biosample"].strip().upper()
    dcc_biosample_clean = dcc_biosample.strip().upper()
    if cand_biosample_clean != dcc_biosample_clean:
        return False, f"REJECTED_CONTRADICTORY_BIOSAMPLE (candidate: '{norm_cand['biosample']}', DCC authority: '{dcc_biosample}')", None, auth_hash

    # Biological perturbation matching (exact identities and baseline rules)
    dcc_treatments = extract_treatment_terms(exp_record)
    cand_treatments = extract_candidate_treatment_terms(norm_cand)
    dcc_mods = extract_genetic_modification_tuples(exp_record)
    cand_mods = extract_candidate_genetic_modification_tuples(norm_cand)

    cand_treat_raw = str(norm_cand.get("treatment_context", "")).strip()
    cand_gm_raw = str(norm_cand.get("genetic_modification_context", "")).strip()

    cand_claims_untreated = (
        not cand_treat_raw or
        "none" in cand_treat_raw.lower() or
        "untreated" in cand_treat_raw.lower()
    )
    cand_claims_unmodified = (
        not cand_gm_raw or
        "none" in cand_gm_raw.lower() or
        "wild-type" in cand_gm_raw.lower() or
        "unmodified" in cand_gm_raw.lower()
    )

    # Contradictory treatments (e.g. dexamethasone vs doxorubicin)
    if dcc_treatments and cand_treatments and dcc_treatments != cand_treatments:
        return False, (
            f"REJECTED_BIOLOGICAL_TREATMENT_OR_GENETIC_MODIFICATION_DISCORDANCE: "
            f"Contradictory biological treatment identity (DCC authority: {sorted(dcc_treatments)}, "
            f"candidate claims: {sorted(cand_treatments)})"
        ), None, auth_hash

    # Contradictory genetic modifications (e.g. CRISPR TP53 vs shRNA RUNX1)
    if dcc_mods and cand_mods and dcc_mods != cand_mods:
        return False, (
            f"REJECTED_BIOLOGICAL_TREATMENT_OR_GENETIC_MODIFICATION_DISCORDANCE: "
            f"Contradictory genetic modification identity (DCC authority: {sorted(dcc_mods)}, "
            f"candidate claims: {sorted(cand_mods)})"
        ), None, auth_hash

    # Direction 1: DCC has perturbations, candidate claims untreated baseline
    if dcc_treatments and cand_claims_untreated:
        return False, (
            f"REJECTED_BIOLOGICAL_TREATMENT_OR_GENETIC_MODIFICATION_DISCORDANCE: "
            f"DCC records treatment {sorted(dcc_treatments)}, but candidate specified untreated baseline"
        ), None, auth_hash

    if dcc_mods and cand_claims_unmodified:
        return False, (
            f"REJECTED_BIOLOGICAL_TREATMENT_OR_GENETIC_MODIFICATION_DISCORDANCE: "
            f"DCC records genetic modification {sorted(dcc_mods)}, but candidate specified wild-type unmodified"
        ), None, auth_hash

    # Direction 2: Candidate claims perturbations, DCC describes untreated baseline
    if cand_treatments and not dcc_treatments:
        return False, (
            f"REJECTED_UNSUPPORTED_CANDIDATE_PERTURBATION_CLAIM: "
            f"Candidate claims treatment {sorted(cand_treatments)}, but authoritative DCC experiment describes untreated baseline"
        ), None, auth_hash

    if cand_mods and not dcc_mods:
        return False, (
            f"REJECTED_UNSUPPORTED_CANDIDATE_PERTURBATION_CLAIM: "
            f"Candidate claims genetic modification {sorted(cand_mods)}, but authoritative DCC experiment describes wild-type unmodified baseline"
        ), None, auth_hash

    # Independent human gene/transcript/TSS annotation verification
    annotator = gene_annotator or get_gene_annotator()
    is_gene_valid, gene_code, gene_info = annotator.validate_gene_and_coordinates(
        symbol=norm_cand["gene_symbol"],
        chromosome=norm_cand["chromosome"],
        interval_start=norm_cand["interval_start"],
        interval_end=norm_cand["interval_end"],
    )
    if not is_gene_valid:
        return False, gene_code, None, auth_hash

    verified_info = {
        "experiment_accession": exp_acc,
        "file_accession": file_acc,
        "tf_target": dcc_tf,
        "biosample": dcc_biosample,
        "assembly": "GRCh38",
        "file_size": dcc_size,
        "md5sum": dcc_md5,
        "canonical_gene": gene_info["symbol"] if gene_info else norm_cand["gene_symbol"],
        "canonical_tss": gene_info["canonical_tss"] if gene_info else None,
        "authority_provenance": authority_store.provenance_source,
    }
    return True, "AUTHENTICATED_SOURCE_IDENTITY", verified_info, auth_hash


# -----------------------------------------------------------------------------
# Candidate Screener Class (Curator Lifecycle & Atomic Ledger)
# -----------------------------------------------------------------------------

class CandidateScreener:
    """
    Curator Candidate Screener for Phase 4C-3 Independent Curators.

    Enforces:
    1. Verifiable ENCODE source authority:
       - Only PRODUCTION_PINNED_SNAPSHOT and PRODUCTION_ENCODE_PORTAL can authenticate.
       - Rejects bare dictionaries, unpinned stores, and synthetic test fixtures fail-closed.
    2. Strict distinction between:
       - AUTHENTICATED_ELIGIBLE_CANDIDATE: both metadata and primary BED bytes verified.
       - PROVISIONAL_METADATA_VERIFIED_CANDIDATE: official metadata verified; BED bytes not yet hashed on disk.
       - UNEXPOSED_BY_REGISTER_CHECK: passed historical register; awaiting authority verification.
       - REJECTED_*: fail-closed rejection.
    3. Collision-safe atomic screening IDs (SCR-XXXX) across concurrent curator instances.
    4. Audited ledger recording dual hashes (candidate metadata SHA-256 and authority response SHA-256).
    """

    LEDGER_COLUMNS = [
        "screening_id",
        "candidate_accession",
        "file_accession",
        "gene_symbol",
        "transcription_factor",
        "biosample_lineage",
        "assembly",
        "genomic_interval",
        "retrieval_source",
        "retrieval_timestamp",
        "curator_id",
        "evaluation_decision",
        "exclusion_reason",
        "source_file_checksum",
        "sha256_candidate_metadata",
        "sha256_authority_response",
    ]

    def __init__(
        self,
        register_path: str | Path | None = None,
        ledger_path: str | Path | None = None,
        authority_store: AuthoritativeDCCStore | Any | None = None,
        authority_catalog: Any | None = None,  # Backwards-compatible argument name
        gene_annotator: GeneAnnotationReference | None = None,
        verify_frozen_hash: bool = True,
        allow_synthetic_fixtures: bool = False,
    ):
        self.register_records = load_authoritative_register(
            register_path=register_path,
            verify_frozen_hash=verify_frozen_hash,
        )
        self.ledger_path = Path(ledger_path) if ledger_path else Path("PHASE4C3_SCREENING_AND_EXCLUSION_LEDGER.csv")
        self.authority_store = authority_store if authority_store is not None else authority_catalog
        self.gene_annotator = gene_annotator or get_gene_annotator()
        self.allow_synthetic_fixtures = allow_synthetic_fixtures

    def screen_candidate(
        self,
        candidate: dict[str, Any],
        curator_id: str = "CURATOR_A",
        write_ledger: bool = True,
        require_authority: bool = True,
        bed_file_path: str | Path | None = None,
        retrieval_source: str = "https://www.encodeproject.org/",
    ) -> ScreeningResult:
        """
        Screens candidate record through fail-closed validation, register purity check,
        authoritative source verification, and independent primary BED byte authentication.
        """
        cand_hash = compute_candidate_hash(candidate)

        # Step 1: Input Validation (Fail-Closed)
        is_valid, val_code, norm_cand = validate_candidate_metadata(candidate)
        if not is_valid or norm_cand is None:
            return self._finalize_result(
                is_admissible=False,
                evaluation_decision="REJECTED_MALFORMED_METADATA",
                reason_code=val_code,
                matched_record_id=None,
                cand_hash=cand_hash,
                auth_hash="NONE_MALFORMED_INPUT",
                cand_summary=candidate if isinstance(candidate, dict) else {},
                curator_id=curator_id,
                retrieval_source=retrieval_source,
                write_ledger=write_ledger,
            )

        # Step 2: Banned Exposure Register Check (Purity & Leakage)
        is_unexposed, reason_code, matched_id = evaluate_candidate_for_leakage(
            norm_cand, self.register_records, verify_frozen_hash=False
        )
        if not is_unexposed:
            return self._finalize_result(
                is_admissible=False,
                evaluation_decision="REJECTED_BANNED_EXPOSURE",
                reason_code=reason_code,
                matched_record_id=matched_id,
                cand_hash=cand_hash,
                auth_hash="NONE_EXCLUDED_BANNED",
                cand_summary=norm_cand,
                curator_id=curator_id,
                retrieval_source=retrieval_source,
                write_ledger=write_ledger,
            )

        # Step 3: Authoritative ENCODE Source Verification
        if require_authority:
            if self.authority_store is None:
                return self._finalize_result(
                    is_admissible=False,
                    evaluation_decision="REJECTED_UNAUTHENTICATED_ACCESSION",
                    reason_code="MISSING_AUTHORITY_CATALOG_CANNOT_AUTHENTICATE",
                    matched_record_id=None,
                    cand_hash=cand_hash,
                    auth_hash="NONE_MISSING_AUTHORITY",
                    cand_summary=norm_cand,
                    curator_id=curator_id,
                    retrieval_source=retrieval_source,
                    write_ledger=write_ledger,
                )

            # Untrusted bare dictionary rejection
            if not isinstance(self.authority_store, AuthoritativeDCCStore):
                return self._finalize_result(
                    is_admissible=False,
                    evaluation_decision="REJECTED_UNTRUSTED_AUTHORITY_SOURCE",
                    reason_code="SUPPLIED_DICTIONARY_LACKS_TRUST_ANCHOR: A supplied dictionary cannot authenticate itself",
                    matched_record_id=None,
                    cand_hash=cand_hash,
                    auth_hash="NONE_UNTRUSTED_DICTIONARY",
                    cand_summary=norm_cand,
                    curator_id=curator_id,
                    retrieval_source=retrieval_source,
                    write_ledger=write_ledger,
                )

            # Synthetic test fixture quarantine: Cannot grant production eligibility
            if self.authority_store.trust_level not in {
                AuthorityTrustLevel.PRODUCTION_PINNED_SNAPSHOT,
                AuthorityTrustLevel.PRODUCTION_ENCODE_PORTAL,
            }:
                if not self.allow_synthetic_fixtures:
                    return self._finalize_result(
                        is_admissible=False,
                        evaluation_decision="REJECTED_SYNTHETIC_TEST_FIXTURE",
                        reason_code=(
                            "SYNTHETIC_FIXTURE_CANNOT_AUTHENTICATE_PRODUCTION_CANDIDATE: "
                            "Test fixtures are synthetic and cannot grant production eligibility"
                        ),
                        matched_record_id=None,
                        cand_hash=cand_hash,
                        auth_hash=self.authority_store.provenance_hash,
                        cand_summary=norm_cand,
                        curator_id=curator_id,
                        retrieval_source=retrieval_source,
                        write_ledger=write_ledger,
                    )

            auth_valid, auth_code, verified_details, auth_hash = verify_authoritative_encode_identity(
                norm_cand=norm_cand,
                authority_store=self.authority_store,
                gene_annotator=self.gene_annotator,
            )
            if not auth_valid:
                return self._finalize_result(
                    is_admissible=False,
                    evaluation_decision="REJECTED_UNAUTHENTICATED_ACCESSION",
                    reason_code=auth_code,
                    matched_record_id=None,
                    cand_hash=cand_hash,
                    auth_hash=auth_hash,
                    cand_summary=norm_cand,
                    curator_id=curator_id,
                    retrieval_source=retrieval_source,
                    write_ledger=write_ledger,
                )

            # Step 4: Primary BED Byte Authentication vs Provisional Metadata-Verified
            full_summary = dict(norm_cand)
            if verified_details:
                full_summary.update(verified_details)

            target_bed_path = bed_file_path or norm_cand.get("bed_file_path")
            if target_bed_path is not None:
                p_bed = Path(target_bed_path)
                if not p_bed.exists():
                    return self._finalize_result(
                        is_admissible=False,
                        evaluation_decision="REJECTED_MISSING_PRIMARY_BED_FILE",
                        reason_code=f"PRIMARY_BED_FILE_NOT_FOUND: {p_bed}",
                        matched_record_id=None,
                        cand_hash=cand_hash,
                        auth_hash=auth_hash,
                        cand_summary=full_summary,
                        curator_id=curator_id,
                        retrieval_source=retrieval_source,
                        write_ledger=write_ledger,
                    )

                bed_bytes = p_bed.read_bytes()
                is_bed_valid, bed_format_err = validate_bed_format(bed_bytes, p_bed.name)
                if not is_bed_valid:
                    return self._finalize_result(
                        is_admissible=False,
                        evaluation_decision="REJECTED_INVALID_BED_FORMAT",
                        reason_code=bed_format_err,
                        matched_record_id=None,
                        cand_hash=cand_hash,
                        auth_hash=auth_hash,
                        cand_summary=full_summary,
                        curator_id=curator_id,
                        retrieval_source=retrieval_source,
                        write_ledger=write_ledger,
                    )

                calc_size = len(bed_bytes)
                calc_md5 = compute_bytes_md5(bed_bytes)
                dcc_size = verified_details.get("file_size") if verified_details else None
                dcc_md5 = verified_details.get("md5sum") if verified_details else None

                if not dcc_size or not dcc_md5 or calc_size != dcc_size or calc_md5.lower() != str(dcc_md5).lower():
                    return self._finalize_result(
                        is_admissible=False,
                        evaluation_decision="REJECTED_PRIMARY_BED_BYTE_CHECKSUM_OR_SIZE_MISMATCH",
                        reason_code=(
                            f"BYTE_VERIFICATION_FAILED: calculated size={calc_size} (expected {dcc_size}), "
                            f"md5={calc_md5} (expected {dcc_md5})"
                        ),
                        matched_record_id=None,
                        cand_hash=cand_hash,
                        auth_hash=auth_hash,
                        cand_summary=full_summary,
                        curator_id=curator_id,
                        retrieval_source=retrieval_source,
                        write_ledger=write_ledger,
                    )

                full_summary["verified_bed_size"] = calc_size
                full_summary["verified_bed_md5"] = calc_md5
                return self._finalize_result(
                    is_admissible=True,
                    evaluation_decision="AUTHENTICATED_ELIGIBLE_CANDIDATE",
                    reason_code="AUTHENTICATED_SOURCE_AND_PRIMARY_BED_BYTES_VERIFIED",
                    matched_record_id=None,
                    cand_hash=cand_hash,
                    auth_hash=auth_hash,
                    cand_summary=full_summary,
                    curator_id=curator_id,
                    retrieval_source=retrieval_source,
                    write_ledger=write_ledger,
                )

            # Metadata verified, but physical primary BED bytes not present on disk
            return self._finalize_result(
                is_admissible=True,
                evaluation_decision="PROVISIONAL_METADATA_VERIFIED_CANDIDATE",
                reason_code="METADATA_VERIFIED_AWAITING_PRIMARY_BED_BYTE_AUTHENTICATION",
                matched_record_id=None,
                cand_hash=cand_hash,
                auth_hash=auth_hash,
                cand_summary=full_summary,
                curator_id=curator_id,
                retrieval_source=retrieval_source,
                write_ledger=write_ledger,
            )

        # Authority check explicitly unperformed
        return self._finalize_result(
            is_admissible=False,
            evaluation_decision="UNEXPOSED_BY_REGISTER_CHECK",
            reason_code="UNEXPOSED_AWAITING_SOURCE_AUTHENTICATION",
            matched_record_id=None,
            cand_hash=cand_hash,
            auth_hash="NONE_AUTHORITY_UNREQUESTED",
            cand_summary=norm_cand,
            curator_id=curator_id,
            retrieval_source=retrieval_source,
            write_ledger=write_ledger,
        )

    def _finalize_result(
        self,
        is_admissible: bool,
        evaluation_decision: str,
        reason_code: str,
        matched_record_id: str | None,
        cand_hash: str,
        auth_hash: str,
        cand_summary: dict[str, Any],
        curator_id: str,
        retrieval_source: str,
        write_ledger: bool,
    ) -> ScreeningResult:
        """Assigns sequential unique screening ID under atomic file lock and writes to ledger."""
        if write_ledger:
            with LedgerFileLock(self.ledger_path):
                max_id = 0
                file_exists = self.ledger_path.exists() and self.ledger_path.stat().st_size > 0
                if file_exists:
                    try:
                        with open(self.ledger_path, "r", encoding="utf-8") as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                sid = row.get("screening_id", "")
                                m = re.match(r"^SCR-(\d+)$", sid)
                                if m:
                                    max_id = max(max_id, int(m.group(1)))
                    except Exception:
                        pass

                next_id = max_id + 1
                screening_id = f"SCR-{next_id:04d}"

                chrom = cand_summary.get("chromosome", "")
                start = cand_summary.get("interval_start", "")
                end = cand_summary.get("interval_end", "")
                genomic_interval = f"{chrom}:{start}-{end}" if chrom and start != "" and end != "" else ""
                checksum = (
                    cand_summary.get("verified_bed_md5") or
                    cand_summary.get("md5sum") or
                    cand_summary.get("authoritative_md5") or
                    ""
                )

                row = {
                    "screening_id": screening_id,
                    "candidate_accession": cand_summary.get("experiment_accession") or "MALFORMED_ACCESSION",
                    "file_accession": cand_summary.get("file_accession") or "",
                    "gene_symbol": cand_summary.get("gene_symbol") or "MALFORMED_GENE",
                    "transcription_factor": cand_summary.get("tf_target") or cand_summary.get("transcription_factor") or "",
                    "biosample_lineage": cand_summary.get("biosample") or cand_summary.get("biosample_lineage") or "",
                    "assembly": cand_summary.get("assembly") or "",
                    "genomic_interval": genomic_interval,
                    "retrieval_source": retrieval_source,
                    "retrieval_timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    "curator_id": curator_id,
                    "evaluation_decision": evaluation_decision,
                    "exclusion_reason": f"{reason_code}{f' (matched: {matched_record_id})' if matched_record_id else ''}",
                    "source_file_checksum": checksum,
                    "sha256_candidate_metadata": cand_hash,
                    "sha256_authority_response": auth_hash,
                }

                with open(self.ledger_path, "a", newline="", encoding="utf-8") as f:
                    writer = csv.DictWriter(f, fieldnames=self.LEDGER_COLUMNS)
                    if not file_exists:
                        writer.writeheader()
                    writer.writerow(row)
                    f.flush()
                    try:
                        os.fsync(f.fileno())
                    except Exception:
                        pass
        else:
            screening_id = "SCR-PREVIEW"

        return ScreeningResult(
            is_admissible=is_admissible,
            evaluation_decision=evaluation_decision,
            reason_code=reason_code,
            matched_record_id=matched_record_id,
            sha256_candidate_metadata=cand_hash,
            sha256_authority_response=auth_hash,
            screening_id=screening_id,
            candidate_summary=cand_summary,
        )
