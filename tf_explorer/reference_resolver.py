"""
Dynamic Authoritative Gene and Transcript Reference Resolver for TF-Explorer.
Phase 4A Generalization Engine.

Resolves arbitrary human genes to verified Ensembl and MANE Select references
without relying on fixed hardcoded gene registries. Fails closed with explicit
categorical status codes when ambiguity or discordance is detected.
"""

from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from enum import Enum
import hashlib
import json
import logging
import os
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
import requests
try:
    from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type
except ImportError:
    def retry(*args, **kwargs):
        def decorator(f):
            return f
        return decorator
    def stop_after_attempt(*args, **kwargs): return None
    def wait_exponential(*args, **kwargs): return None
    def retry_if_exception_type(*args, **kwargs): return None

logger = logging.getLogger(__name__)

ENSEMBL_REST_URL = "https://rest.ensembl.org"
DEFAULT_TIMEOUT = 10


class ReferenceResolverStatus(str, Enum):
    """Categorical response states for reference resolution."""
    REFERENCE_VERIFIED = "REFERENCE_VERIFIED"
    GENE_AMBIGUOUS = "GENE_AMBIGUOUS"
    TRANSCRIPT_NOT_FOUND = "TRANSCRIPT_NOT_FOUND"
    TRANSCRIPT_VERSION_MISMATCH = "TRANSCRIPT_VERSION_MISMATCH"
    ASSEMBLY_MISMATCH = "ASSEMBLY_MISMATCH"
    REFERENCE_UNAVAILABLE = "REFERENCE_UNAVAILABLE"
    REFERENCE_METADATA_INCOMPLETE = "REFERENCE_METADATA_INCOMPLETE"
    REFERENCE_SOURCE_CONFLICT = "REFERENCE_SOURCE_CONFLICT"
    MANE_UNVERIFIED = "MANE_UNVERIFIED"
    MANE_SELECT_VERIFIED = "MANE_SELECT_VERIFIED"
    MANE_PLUS_CLINICAL_VERIFIED = "MANE_PLUS_CLINICAL_VERIFIED"
    MANE_VERSION_MISMATCH = "MANE_VERSION_MISMATCH"


class EnsemblResolverRateLimitError(Exception):
    """Raised when Ensembl API returns HTTP 429 rate limit."""
    pass


class EnsemblResolverServerError(Exception):
    """Raised when Ensembl API returns HTTP 5xx server error."""
    pass


@dataclass
class TranscriptReference:
    """Detailed model of a single transcript reference."""
    transcript_id: str
    transcript_version: str
    display_name: str
    biotype: str
    is_canonical: bool
    is_mane_select: bool
    mane_select_id: Optional[str]
    chromosome: str
    strand: str
    start: int
    end: int
    tss: int
    length_bp: int
    exons_count: Optional[int] = None
    is_mane_plus_clinical: bool = False
    mane_status: str = "MANE_UNVERIFIED"
    mane_release: str = "MANE Release 1.3"

    @property
    def version(self) -> str:
        return self.transcript_version

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class GeneReference:
    """Complete authoritative reference record for an arbitrary resolved gene."""
    gene_symbol: str
    ensembl_gene_id: str
    ncbi_gene_id: Optional[str]
    species: str
    genome_assembly: str
    chromosome: str
    strand: str
    gene_start: int
    gene_end: int
    description: Optional[str]
    transcripts: List[TranscriptReference]
    canonical_transcript: Optional[TranscriptReference]
    mane_select_transcript: Optional[TranscriptReference]
    annotation_release: str
    source_url: str
    retrieval_timestamp: str
    status: ReferenceResolverStatus
    status_notes: List[str] = field(default_factory=list)
    mane_plus_clinical_transcript: Optional[TranscriptReference] = None
    mane_status: str = "MANE_UNVERIFIED"
    mane_release: str = "MANE Release 1.3"

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["status"] = self.status.value
        d["transcripts"] = [t.to_dict() for t in self.transcripts]
        d["canonical_transcript"] = self.canonical_transcript.to_dict() if self.canonical_transcript else None
        d["mane_select_transcript"] = self.mane_select_transcript.to_dict() if self.mane_select_transcript else None
        d["mane_plus_clinical_transcript"] = self.mane_plus_clinical_transcript.to_dict() if self.mane_plus_clinical_transcript else None
        return d


@dataclass
class ResolutionResult:
    """Outcome of attempting to resolve a gene and transcript reference."""
    status: ReferenceResolverStatus
    gene_reference: Optional[GeneReference] = None
    ambiguous_candidates: List[Dict[str, Any]] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def is_verified(self) -> bool:
        return self.status == ReferenceResolverStatus.REFERENCE_VERIFIED and self.gene_reference is not None

    @property
    def canonical_transcript(self) -> Optional[TranscriptReference]:
        return self.gene_reference.canonical_transcript if self.gene_reference else None

    @property
    def canonical_tss(self) -> Optional[int]:
        if self.canonical_transcript:
            return self.canonical_transcript.tss
        return None

    @property
    def gene_symbol(self) -> Optional[str]:
        return self.gene_reference.gene_symbol if self.gene_reference else None

    @property
    def gene_start(self) -> Optional[int]:
        return self.gene_reference.gene_start if self.gene_reference else None

    @property
    def gene_end(self) -> Optional[int]:
        return self.gene_reference.gene_end if self.gene_reference else None

    @property
    def strand(self) -> Optional[str]:
        return self.gene_reference.strand if self.gene_reference else None

    @property
    def chromosome(self) -> Optional[str]:
        return self.gene_reference.chromosome if self.gene_reference else None

    @property
    def transcripts(self) -> List[TranscriptReference]:
        return self.gene_reference.transcripts if self.gene_reference else []

    @property
    def ensembl_gene_id(self) -> Optional[str]:
        return self.gene_reference.ensembl_gene_id if self.gene_reference else None

    @property
    def ncbi_gene_id(self) -> Optional[str]:
        return self.gene_reference.ncbi_gene_id if self.gene_reference else None

    @property
    def error_message(self) -> str:
        return "; ".join(self.errors) if self.errors else ""

    def get_promoter_interval_1based(self, upstream: int = 2000, downstream: int = 500) -> Tuple[int, int]:
        """Calculates 1-based inclusive promoter interval for canonical transcript."""
        if not self.is_verified or not self.canonical_transcript:
            raise ValueError("Cannot calculate promoter interval on unverified or missing reference")
        tss = self.canonical_transcript.tss
        strand = self.strand
        if strand == "+":
            return (tss - upstream, tss + downstream)
        elif strand == "-":
            return (tss - downstream, tss + upstream)
        else:
            raise ValueError(f"Unknown strand: {strand}")

    def get_promoter_interval_bed0based(self, upstream: int = 2000, downstream: int = 500) -> Tuple[int, int]:
        """Calculates 0-based half-open (BED) interval for canonical transcript."""
        p1 = self.get_promoter_interval_1based(upstream, downstream)
        return (p1[0] - 1, p1[1])

    def to_dict(self) -> Dict[str, Any]:
        return {
            "status": self.status.value,
            "is_verified": self.is_verified,
            "gene_reference": self.gene_reference.to_dict() if self.gene_reference else None,
            "ambiguous_candidates": self.ambiguous_candidates,
            "errors": self.errors,
            "warnings": self.warnings
        }


class DynamicReferenceResolver:
    """
    Generalizable, reference-aware gene and transcript resolver.
    Connects gene identity, Ensembl annotations, MANE Select transcripts,
    and assembly strand orientation without hardcoded gene lists.
    """

    REGEX_ENSEMBL_GENE = re.compile(r"^ENSG[0-9]{11}(\.[0-9]+)?$")
    REGEX_ENSEMBL_TRANSCRIPT = re.compile(r"^ENST[0-9]{11}(\.[0-9]+)?$")

    def __init__(
        self,
        cache_dir: Optional[str] = None,
        offline_snapshots: Optional[Dict[str, Any]] = None,
        snapshot_file: Optional[Union[str, Path]] = None,
        mane_file: Optional[Union[str, Path]] = None,
        timeout: int = DEFAULT_TIMEOUT
    ):
        self.cache_dir = str(cache_dir) if cache_dir else os.path.join("reference_data", "cache")
        self.offline_snapshots = dict(offline_snapshots) if offline_snapshots else {}
        self.timeout = timeout
        os.makedirs(self.cache_dir, exist_ok=True)

        self.mane_select_index: Dict[str, Any] = {}
        if mane_file:
            m_path = Path(mane_file)
            if m_path.exists():
                try:
                    with open(m_path, "r", encoding="utf-8") as f:
                        self.mane_select_index = json.load(f)
                except Exception as e_m:
                    logger.warning(f"Could not load MANE index from {m_path}: {e_m}")
        elif offline_snapshots is None or snapshot_file is not None:
            default_mane = Path(__file__).resolve().parent.parent / "reference_data" / "MANE_v1.3_select.json"
            if default_mane.exists():
                try:
                    with open(default_mane, "r", encoding="utf-8") as f:
                        self.mane_select_index = json.load(f)
                except Exception as e_m:
                    logger.warning(f"Could not load default MANE index: {e_m}")

        # Seed offline snapshots from snapshot_file or default reference_data/gene_snapshots.json
        if snapshot_file:
            snap_path = str(snapshot_file)
            if os.path.exists(snap_path):
                try:
                    with open(snap_path, "r", encoding="utf-8") as f:
                        snaps = json.load(f)
                        for k, v in snaps.items():
                            self.offline_snapshots.setdefault(k, v)
                except Exception as e_sn:
                    logger.warning(f"Could not load gene snapshots from {snap_path}: {e_sn}")
        elif offline_snapshots is None:
            default_snap = os.path.join("reference_data", "gene_snapshots.json")
            if os.path.exists(default_snap):
                try:
                    with open(default_snap, "r", encoding="utf-8") as f:
                        snaps = json.load(f)
                        for k, v in snaps.items():
                            self.offline_snapshots.setdefault(k, v)
                except Exception as e_sn:
                    logger.warning(f"Could not load default gene snapshots: {e_sn}")

    @retry(
        reraise=True,
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=1, max=6),
        retry=retry_if_exception_type((
            requests.exceptions.ConnectionError,
            requests.exceptions.Timeout,
            EnsemblResolverRateLimitError,
            EnsemblResolverServerError
        ))
    )
    def _http_get_json(self, url: str) -> requests.Response:
        """Executes GET request with proper headers, user-agent, and error handling."""
        headers = {
            "Content-Type": "application/json",
            "User-Agent": "TF-Explorer/4.0.0 (Research Software; Ensembl-REST)"
        }
        resp = requests.get(url, headers=headers, timeout=self.timeout)
        if resp.status_code == 429:
            retry_after = resp.headers.get("Retry-After", "1.0")
            try:
                time.sleep(float(retry_after))
            except ValueError:
                time.sleep(1.0)
            raise EnsemblResolverRateLimitError(f"Ensembl rate limit reached: {url}")
        if resp.status_code >= 500:
            raise EnsemblResolverServerError(f"Ensembl server error HTTP {resp.status_code}: {url}")
        return resp

    def resolve(
        self,
        query: str,
        assembly: str = "GRCh38",
        requested_transcript_id: Optional[str] = None,
        requested_transcript_version: Optional[str] = None,
        transcript_id: Optional[str] = None,
        transcript_version: Optional[str] = None,
        require_mane: bool = False
    ) -> ResolutionResult:
        """
        Main entry point. Resolves either an Ensembl gene ID (ENSG...) or symbol.
        """
        requested_transcript_id = requested_transcript_id or transcript_id
        requested_transcript_version = requested_transcript_version or transcript_version
        if requested_transcript_id and "." in requested_transcript_id and not requested_transcript_version:
            parts = requested_transcript_id.split(".", 1)
            requested_transcript_id = parts[0]
            requested_transcript_version = parts[1]
        clean_query = query.strip()
        if not clean_query:
            return ResolutionResult(
                status=ReferenceResolverStatus.REFERENCE_METADATA_INCOMPLETE,
                errors=["Gene query string is empty."]
            )

        # 1. Check assembly support (currently restricted to GRCh38)
        if assembly.upper() not in ["GRCH38", "HG38"]:
            return ResolutionResult(
                status=ReferenceResolverStatus.ASSEMBLY_MISMATCH,
                errors=[f"Requested assembly '{assembly}' is not supported. Only GRCh38 is currently supported."]
            )

        # 2. Check offline snapshot dictionary first (for hermetic test execution)
        norm_key = clean_query.upper()
        if norm_key in self.offline_snapshots:
            snapshot_data = self.offline_snapshots[norm_key]
            return self._parse_ensembl_payload(
                snapshot_data,
                clean_query,
                assembly="GRCh38",
                requested_transcript_id=requested_transcript_id,
                requested_transcript_version=requested_transcript_version,
                source_url="offline_snapshot",
                require_mane=require_mane
            )

        # 3. Check persistent disk cache
        cache_key = hashlib.sha256(f"{norm_key}_{assembly}".encode("utf-8")).hexdigest()
        cache_file = os.path.join(self.cache_dir, f"{cache_key}.json")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached = json.load(f)
                    return self._parse_ensembl_payload(
                        cached["payload"],
                        clean_query,
                        assembly="GRCh38",
                        requested_transcript_id=requested_transcript_id,
                        requested_transcript_version=requested_transcript_version,
                        source_url=cached.get("source_url", "cached_response"),
                        require_mane=require_mane
                    )
            except Exception as e_cache:
                logger.warning(f"Failed reading cache for {clean_query}: {e_cache}")

        # 4. Resolve via Ensembl REST API with expand=1;mane=1
        if self.REGEX_ENSEMBL_GENE.match(clean_query):
            # Query by Ensembl ID directly
            url = f"{ENSEMBL_REST_URL}/lookup/id/{clean_query}?expand=1;mane=1"
            return self._fetch_and_parse(url, clean_query, assembly, requested_transcript_id, requested_transcript_version, cache_file, require_mane=require_mane)
        else:
            # Query by Symbol
            # First check symbol exactness
            url = f"{ENSEMBL_REST_URL}/lookup/symbol/homo_sapiens/{clean_query}?expand=1;mane=1"
            return self._fetch_and_parse(url, clean_query, assembly, requested_transcript_id, requested_transcript_version, cache_file, require_mane=require_mane)

    def _fetch_and_parse(
        self,
        url: str,
        query: str,
        assembly: str,
        requested_transcript_id: Optional[str],
        requested_transcript_version: Optional[str],
        cache_file: str,
        require_mane: bool = False
    ) -> ResolutionResult:
        try:
            resp = self._http_get_json(url)
            if resp.status_code == 404:
                # Symbol not found or ambiguous
                return self._handle_not_found(query)

            if not resp.ok:
                return ResolutionResult(
                    status=ReferenceResolverStatus.REFERENCE_UNAVAILABLE,
                    errors=[f"Ensembl REST API returned HTTP {resp.status_code}: {resp.text}"]
                )

            payload = resp.json()
            if not isinstance(payload, dict) or "id" not in payload:
                return ResolutionResult(
                    status=ReferenceResolverStatus.REFERENCE_METADATA_INCOMPLETE,
                    errors=[f"Incomplete payload returned by Ensembl for query '{query}'."]
                )

            # Save clean payload to cache
            try:
                with open(cache_file, "w", encoding="utf-8") as fh:
                    json.dump({"payload": payload, "source_url": url, "timestamp": datetime.now(timezone.utc).isoformat()}, fh)
            except Exception as e_w:
                logger.warning(f"Could not write cache file {cache_file}: {e_w}")

            return self._parse_ensembl_payload(
                payload,
                query,
                assembly,
                requested_transcript_id,
                requested_transcript_version,
                source_url=url,
                require_mane=require_mane
            )

        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e_net:
            logger.error(f"Network error resolving {query}: {e_net}")
            return ResolutionResult(
                status=ReferenceResolverStatus.REFERENCE_UNAVAILABLE,
                errors=[f"Network communication failed while contacting Ensembl REST API: {str(e_net)}"]
            )
        except Exception as e:
            logger.error(f"Unexpected error resolving {query}: {e}")
            return ResolutionResult(
                status=ReferenceResolverStatus.REFERENCE_UNAVAILABLE,
                errors=[f"Encountered unexpected resolution error: {str(e)}"]
            )

    def _handle_not_found(self, query: str) -> ResolutionResult:
        """Attempts xrefs symbol lookup to distinguish unknown from ambiguous symbols."""
        xrefs_url = f"{ENSEMBL_REST_URL}/xrefs/symbol/homo_sapiens/{query}"
        try:
            r = self._http_get_json(xrefs_url)
            if r.ok:
                candidates = r.json()
                if isinstance(candidates, list) and len(candidates) > 1:
                    gene_hits = [c for c in candidates if c.get("type") == "gene"]
                    if len(gene_hits) > 1:
                        return ResolutionResult(
                            status=ReferenceResolverStatus.GENE_AMBIGUOUS,
                            ambiguous_candidates=gene_hits,
                            errors=[f"Gene symbol '{query}' is ambiguous: matches {len(gene_hits)} distinct Ensembl gene records."]
                        )
        except Exception:
            pass

        return ResolutionResult(
            status=ReferenceResolverStatus.REFERENCE_UNAVAILABLE,
            errors=[f"Gene '{query}' was not found in Ensembl GRCh38 release annotations."]
        )

    def _parse_ensembl_payload(
        self,
        payload: Dict[str, Any],
        query: str,
        assembly: str,
        requested_transcript_id: Optional[str],
        requested_transcript_version: Optional[str],
        source_url: str,
        require_mane: bool = False
    ) -> ResolutionResult:
        """
        Parses authoritative Ensembl payload into strict GeneReference and TranscriptReference objects.
        """
        errors: List[str] = []
        warnings: List[str] = []

        gene_symbol = payload.get("display_name", "").strip() or query
        ensembl_gene_id = payload.get("id", "").strip()
        payload_assembly = payload.get("assembly_name", "GRCh38")

        if payload_assembly.upper() != "GRCH38":
            return ResolutionResult(
                status=ReferenceResolverStatus.ASSEMBLY_MISMATCH,
                errors=[f"Payload assembly is '{payload_assembly}', expected 'GRCh38'."]
            )

        seq_region = str(payload.get("seq_region_name", ""))
        chromosome = f"chr{seq_region}" if not seq_region.startswith("chr") else seq_region
        raw_strand = payload.get("strand", 1)
        strand = "+" if raw_strand == 1 or raw_strand == "+" else "-"
        gene_start = int(payload.get("start", 0))
        gene_end = int(payload.get("end", 0))

        # Parse transcripts
        raw_transcripts = payload.get("Transcript", [])
        if not raw_transcripts:
            return ResolutionResult(
                status=ReferenceResolverStatus.REFERENCE_METADATA_INCOMPLETE,
                errors=[f"No transcripts found for gene '{gene_symbol}' ({ensembl_gene_id})."]
            )

        transcripts: List[TranscriptReference] = []
        canonical_tx: Optional[TranscriptReference] = None
        mane_tx: Optional[TranscriptReference] = None
        mane_plus_clinical_tx: Optional[TranscriptReference] = None

        canonical_id_full = payload.get("canonical_transcript", "")
        # canonical_transcript can be "ENST00000257904.11"
        canonical_base_id = canonical_id_full.split(".")[0] if canonical_id_full else ""

        for raw_tx in raw_transcripts:
            t_id = raw_tx.get("id", "")
            t_ver = str(raw_tx.get("version", ""))
            t_display = raw_tx.get("display_name", t_id)
            t_biotype = raw_tx.get("biotype", "unknown")
            t_strand_raw = raw_tx.get("strand", raw_strand)
            t_strand = "+" if t_strand_raw == 1 or t_strand_raw == "+" else "-"
            t_start = int(raw_tx.get("start", 0))
            t_end = int(raw_tx.get("end", 0))

            # Biological TSS Rule:
            # For + strand: 5' TSS is the lower coordinate (start)
            # For - strand: 5' TSS is the higher coordinate (end)
            tss = t_start if t_strand == "+" else t_end

            is_canon = (raw_tx.get("is_canonical") == 1 or t_id == canonical_base_id)

            # Detect MANE status - strictly evidence-backed, NO canonical inference
            mane_raw = raw_tx.get("MANE") or raw_tx.get("MANE_Select") or raw_tx.get("mane_select") or raw_tx.get("mane")
            mane_id = None
            is_mane = False
            is_mane_plus_clinical = False
            tx_mane_status = "MANE_UNVERIFIED"

            # Check official pinned MANE Select Release 1.3 index first
            mane_pinned = self.mane_select_index.get(gene_symbol) or self.mane_select_index.get(query)
            if mane_pinned:
                pinned_ensembl_full = str(mane_pinned.get("Ensembl_nuc", "")).strip()
                pinned_tx_base = pinned_ensembl_full.split(".")[0]
                pinned_tx_ver = pinned_ensembl_full.split(".")[1] if "." in pinned_ensembl_full else ""
                
                if t_id == pinned_tx_base:
                    # Base ID matches! Check transcript version concordance
                    if t_ver and pinned_tx_ver:
                        if t_ver == pinned_tx_ver:
                            is_mane = True
                            tx_mane_status = "MANE_SELECT_VERIFIED"
                            mane_id = mane_pinned.get("RefSeq_nuc")
                        else:
                            # Outdated / divergent transcript version!
                            is_mane = False
                            tx_mane_status = "MANE_VERSION_MISMATCH"
                            mane_id = None
                    else:
                        is_mane = True
                        tx_mane_status = "MANE_SELECT_VERIFIED"
                        mane_id = mane_pinned.get("RefSeq_nuc")

            # Check Ensembl MANE payload structures (distinguishing MANE_Select vs MANE_Plus_Clinical)
            if not is_mane and tx_mane_status != "MANE_VERSION_MISMATCH" and mane_raw:
                if isinstance(mane_raw, list):
                    for item in mane_raw:
                        if isinstance(item, dict):
                            m_type = str(item.get("type", "")).strip().lower()
                            if "plus_clinical" in m_type or "plus" in m_type:
                                is_mane_plus_clinical = True
                                tx_mane_status = "MANE_PLUS_CLINICAL_VERIFIED"
                                mane_id = item.get("refseq_match") or item.get("refseq") or item.get("id")
                            elif "select" in m_type:
                                is_mane = True
                                tx_mane_status = "MANE_SELECT_VERIFIED"
                                mane_id = item.get("refseq_match") or item.get("refseq") or item.get("id")
                                break
                elif isinstance(mane_raw, dict):
                    m_type = str(mane_raw.get("type", "")).strip().lower()
                    if "plus_clinical" in m_type or "plus" in m_type:
                        is_mane_plus_clinical = True
                        tx_mane_status = "MANE_PLUS_CLINICAL_VERIFIED"
                        mane_id = mane_raw.get("refseq_match") or mane_raw.get("refseq") or mane_raw.get("id")
                    elif "select" in m_type:
                        is_mane = True
                        tx_mane_status = "MANE_SELECT_VERIFIED"
                        mane_id = mane_raw.get("refseq_match") or mane_raw.get("refseq") or mane_raw.get("id")
                elif isinstance(mane_raw, str) and mane_raw.strip():
                    s_clean = mane_raw.strip().lower()
                    if "plus_clinical" in s_clean or "plus" in s_clean:
                        is_mane_plus_clinical = True
                        tx_mane_status = "MANE_PLUS_CLINICAL_VERIFIED"
                        mane_id = mane_raw.strip()
                    elif "select" in s_clean:
                        is_mane = True
                        tx_mane_status = "MANE_SELECT_VERIFIED"
                        mane_id = mane_raw.strip()

            if not is_mane and not is_mane_plus_clinical and tx_mane_status != "MANE_VERSION_MISMATCH":
                tags = raw_tx.get("tag", []) if isinstance(raw_tx.get("tag"), list) else []
                flags = raw_tx.get("flags", []) if isinstance(raw_tx.get("flags"), list) else []
                combined_flags = [str(x).lower() for x in tags + flags]
                if any("mane_plus_clinical" in f or "plus_clinical" in f for f in combined_flags):
                    is_mane_plus_clinical = True
                    tx_mane_status = "MANE_PLUS_CLINICAL_VERIFIED"
                elif any("mane_select" in f for f in combined_flags):
                    is_mane = True
                    tx_mane_status = "MANE_SELECT_VERIFIED"

            tx_ref = TranscriptReference(
                transcript_id=t_id,
                transcript_version=t_ver,
                display_name=t_display,
                biotype=t_biotype,
                is_canonical=is_canon,
                is_mane_select=is_mane,
                mane_select_id=mane_id if isinstance(mane_id, str) else None,
                chromosome=chromosome,
                strand=t_strand,
                start=t_start,
                end=t_end,
                tss=tss,
                length_bp=t_end - t_start,
                exons_count=len(raw_tx.get("Exon", [])) if "Exon" in raw_tx else None,
                is_mane_plus_clinical=is_mane_plus_clinical,
                mane_status=tx_mane_status,
                mane_release="MANE Release 1.3"
            )
            transcripts.append(tx_ref)

            if is_canon and canonical_tx is None:
                canonical_tx = tx_ref
            if is_mane and mane_tx is None:
                mane_tx = tx_ref
            if is_mane_plus_clinical and mane_plus_clinical_tx is None:
                mane_plus_clinical_tx = tx_ref

        # Fallback for canonical if not explicitly flagged: sort protein coding then longest
        if canonical_tx is None and transcripts:
            transcripts_sorted = sorted(transcripts, key=lambda x: (x.biotype == "protein_coding", x.length_bp), reverse=True)
            canonical_tx = transcripts_sorted[0]
            canonical_tx.is_canonical = True
            warnings.append(f"Canonical transcript flag missing in source; inferred '{canonical_tx.transcript_id}' as canonical.")

        # If requested transcript specified, verify existence and version
        if requested_transcript_id:
            # Strip version if provided in ID (e.g. ENST00000257904.11)
            if "." in requested_transcript_id:
                parts = requested_transcript_id.split(".", 1)
                target_base_id = parts[0]
                if not requested_transcript_version:
                    requested_transcript_version = parts[1]
            else:
                target_base_id = requested_transcript_id
            matched_tx = next((t for t in transcripts if t.transcript_id == target_base_id), None)
            if not matched_tx:
                return ResolutionResult(
                    status=ReferenceResolverStatus.TRANSCRIPT_NOT_FOUND,
                    errors=[f"Requested transcript '{requested_transcript_id}' was not found in gene '{gene_symbol}' transcripts."]
                )
            if requested_transcript_version and matched_tx.transcript_version != requested_transcript_version:
                return ResolutionResult(
                    status=ReferenceResolverStatus.TRANSCRIPT_VERSION_MISMATCH,
                    errors=[
                        f"Transcript version mismatch for '{target_base_id}': Requested version '{requested_transcript_version}', "
                        f"but Ensembl GRCh38 authoritative version is '{matched_tx.transcript_version}'."
                    ]
                )
            canonical_tx = matched_tx

        # Retrieve or infer NCBI Entrez Gene ID if present in payload or metadata
        ncbi_gene_id = None
        if "entrezgene" in payload:
            ncbi_gene_id = str(payload["entrezgene"])
        elif "ncbi_gene_id" in payload:
            ncbi_gene_id = str(payload["ncbi_gene_id"])

        # Determine overall gene-level MANE status
        if mane_tx is not None:
            overall_mane_status = "MANE_SELECT_VERIFIED"
        elif any(t.mane_status == "MANE_VERSION_MISMATCH" for t in transcripts):
            overall_mane_status = "MANE_VERSION_MISMATCH"
        elif mane_plus_clinical_tx is not None:
            overall_mane_status = "MANE_PLUS_CLINICAL_VERIFIED"
        else:
            overall_mane_status = "MANE_UNVERIFIED"

        if require_mane:
            if overall_mane_status == "MANE_SELECT_VERIFIED":
                pass
            elif overall_mane_status == "MANE_VERSION_MISMATCH":
                mismatched_t = next(t for t in transcripts if t.mane_status == "MANE_VERSION_MISMATCH")
                return ResolutionResult(
                    status=ReferenceResolverStatus.MANE_VERSION_MISMATCH,
                    gene_reference=None,
                    errors=[
                        f"MANE version mismatch for gene '{gene_symbol}': transcript '{mismatched_t.transcript_id}' "
                        f"version '{mismatched_t.transcript_version}' does not match official MANE Release 1.3 specification."
                    ]
                )
            elif overall_mane_status == "MANE_PLUS_CLINICAL_VERIFIED":
                return ResolutionResult(
                    status=ReferenceResolverStatus.MANE_PLUS_CLINICAL_VERIFIED,
                    gene_reference=None,
                    errors=[
                        f"Gene '{gene_symbol}' has MANE Plus Clinical annotation, which is distinct from MANE Select."
                    ]
                )
            else:
                return ResolutionResult(
                    status=ReferenceResolverStatus.MANE_UNVERIFIED,
                    gene_reference=None,
                    errors=[f"MANE Select transcript unverified for gene '{gene_symbol}': No explicit MANE Select annotation present in primary source."]
                )

        status_notes = ["Successfully resolved against authoritative Ensembl GRCh38 metadata."]
        status_notes.append(f"MANE Status: {overall_mane_status} (MANE Release 1.3)")
        if overall_mane_status == "MANE_VERSION_MISMATCH":
            status_notes.append("Historical transcript version detected; does not match current MANE Release 1.3.")
        elif overall_mane_status == "MANE_PLUS_CLINICAL_VERIFIED":
            status_notes.append("Designated as MANE Plus Clinical (distinct from MANE Select).")
        elif overall_mane_status == "MANE_UNVERIFIED":
            status_notes.append("MANE Select transcript unverified: No explicit MANE Select record present in source annotation.")

        gene_ref = GeneReference(
            gene_symbol=gene_symbol,
            ensembl_gene_id=ensembl_gene_id,
            ncbi_gene_id=ncbi_gene_id,
            species="Homo sapiens",
            genome_assembly="GRCh38",
            chromosome=chromosome,
            strand=strand,
            gene_start=gene_start,
            gene_end=gene_end,
            description=payload.get("description"),
            transcripts=transcripts,
            canonical_transcript=canonical_tx,
            mane_select_transcript=mane_tx,
            annotation_release="Ensembl v112 / GRCh38.p14",
            source_url=source_url,
            retrieval_timestamp=datetime.now(timezone.utc).isoformat(),
            status=ReferenceResolverStatus.REFERENCE_VERIFIED,
            status_notes=status_notes,
            mane_plus_clinical_transcript=mane_plus_clinical_tx,
            mane_status=overall_mane_status,
            mane_release="MANE Release 1.3"
        )

        return ResolutionResult(
            status=ReferenceResolverStatus.REFERENCE_VERIFIED,
            gene_reference=gene_ref,
            warnings=warnings
        )

    @staticmethod
    def convert_to_bed_coordinates(start_1based: int, end_1based: int) -> Tuple[int, int]:
        """
        Converts 1-based fully-closed genomic coordinates [start, end]
        to standard BED 0-based half-open coordinates [chromStart, chromEnd).
        Rule: chromStart = start - 1, chromEnd = end.
        """
        if start_1based > end_1based:
            raise ValueError(f"Invalid 1-based interval: start={start_1based} > end={end_1based}.")
        return start_1based - 1, end_1based

    @staticmethod
    def calculate_promoter_coordinates(
        tss: int,
        strand: str,
        upstream: int,
        downstream: int,
        coordinate_system: str = "1-based"
    ) -> Tuple[int, int]:
        """
        Calculates exact promoter interval for a given 1-based TSS.
        
        Strand logic:
        For '+' strand:
          upstream extends towards lower genomic coordinates.
          downstream extends towards higher genomic coordinates.
          Interval (1-based): [tss - upstream, tss + downstream]
        
        For '-' strand:
          upstream extends towards higher genomic coordinates.
          downstream extends towards lower genomic coordinates.
          Interval (1-based): [tss - downstream, tss + upstream]
        """
        if strand not in ["+", "-"]:
            raise ValueError(f"Invalid strand '{strand}'. Must be '+' or '-'.")
        if upstream < 0 or downstream < 0:
            raise ValueError(f"Upstream ({upstream}) and downstream ({downstream}) parameters must be non-negative.")

        if strand == "+":
            p_start = tss - upstream
            p_end = tss + downstream
        else:
            p_start = tss - downstream
            p_end = tss + upstream

        if coordinate_system == "0-based":
            return p_start - 1, p_end
        return p_start, p_end
