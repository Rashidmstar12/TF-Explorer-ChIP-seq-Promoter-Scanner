"""
Evidence Integrity and Dataset Provenance Verification Engine for TF-Explorer.
Phase 3B / Phase 3B-R Enhanced Implementation.

Provides structured, automated verification of ChIP-seq and regulatory genomics datasets
against authoritative ENCODE schemas and standards.
"""

import os
import re
import json
import time
import hashlib
import logging
from enum import Enum
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
import requests

from .provenance import (
    validate_accession_format,
    VERIFIED_METADATA_SNAPSHOT,
    REGEX_ENCODE_FILE,
    REGEX_ENCODE_EXP
)

logger = logging.getLogger(__name__)

ENCODE_BASE_URL = "https://www.encodeproject.org"
CACHE_SCHEMA_VERSION = "2.0"

# Canonical Aliases
TF_ALIASES = {
    "CREB": "CREB1",
    "P53": "TP53",
    "C-MYC": "MYC",
    "ER": "ESR1",
    "ER-ALPHA": "ESR1",
    "ER_ALPHA": "ESR1",
    "ERALPHA": "ESR1",
    "E2F-1": "E2F1"
}

ASSEMBLY_ALIASES = {
    "HG38": "GRCh38",
    "GRCH38": "GRCh38",
    "HG19": "GRCh37",
    "GRCH37": "GRCh37",
    "MM10": "mm10",
    "GRCM38": "mm10",
    "MM39": "GRCm39",
    "GRCM39": "GRCm39"
}


class ValidationStatusCode(str, Enum):
    VALIDATED = "VALIDATED"
    FILE_NOT_FOUND = "FILE_NOT_FOUND"
    EXPERIMENT_NOT_FOUND = "EXPERIMENT_NOT_FOUND"
    FILE_EXPERIMENT_MISMATCH = "FILE_EXPERIMENT_MISMATCH"
    TF_TARGET_MISMATCH = "TF_TARGET_MISMATCH"
    BIOSAMPLE_MISMATCH = "BIOSAMPLE_MISMATCH"
    ASSEMBLY_MISMATCH = "ASSEMBLY_MISMATCH"
    ASSAY_TYPE_MISMATCH = "ASSAY_TYPE_MISMATCH"
    OUTPUT_TYPE_MISMATCH = "OUTPUT_TYPE_MISMATCH"
    STATUS_DEPRECATED_OR_REVOKED = "STATUS_DEPRECATED_OR_REVOKED"
    CHECKSUM_MISMATCH = "CHECKSUM_MISMATCH"
    INVALID_ACCESSION_FORMAT = "INVALID_ACCESSION_FORMAT"
    FILE_UNAVAILABLE = "FILE_UNAVAILABLE"
    METADATA_UNVERIFIED = "METADATA_UNVERIFIED"
    RATE_LIMIT_EXCEEDED = "RATE_LIMIT_EXCEEDED"
    SERVER_ERROR = "SERVER_ERROR"
    CORRUPTED_COORDINATES = "CORRUPTED_COORDINATES"


class ChecksumStatusCode(str, Enum):
    LOCAL_FILE_CHECKSUM_VERIFIED = "LOCAL_FILE_CHECKSUM_VERIFIED"
    REMOTE_CHECKSUM_RETRIEVED = "REMOTE_CHECKSUM_RETRIEVED"
    CHECKSUM_MISMATCH = "CHECKSUM_MISMATCH"
    FILE_NOT_DOWNLOADED = "FILE_NOT_DOWNLOADED"


class ExperimentalContextStatus(str, Enum):
    """Categorical evaluation of biological experimental context compatibility."""
    CONTEXT_MATCH = "CONTEXT_MATCH"
    CONTEXT_MISMATCH = "CONTEXT_MISMATCH"
    CONTEXT_ACCEPTED_WITH_LIMITATIONS = "CONTEXT_ACCEPTED_WITH_LIMITATIONS"
    CONTEXT_UNSPECIFIED = "CONTEXT_UNSPECIFIED"
    CONTEXT_METADATA_INCOMPLETE = "CONTEXT_METADATA_INCOMPLETE"


@dataclass
class ValidationResult:
    """Structured machine-readable evidence verification record."""
    status: ValidationStatusCode
    is_valid: bool
    failure_code: Optional[str] = None
    checksum_status: ChecksumStatusCode = ChecksumStatusCode.FILE_NOT_DOWNLOADED
    metadata_source: str = "ENCODE REST API v2"
    retrieval_date: str = field(default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%d"))
    accessions: Dict[str, str] = field(default_factory=dict)
    message: str = ""
    evidence: Dict[str, Any] = field(default_factory=dict)

    # Biological Experimental Context Verification (Phase 3C-A)
    context_status: ExperimentalContextStatus = ExperimentalContextStatus.CONTEXT_UNSPECIFIED
    context_compatible: bool = True
    context_accepted: bool = False
    context_warnings: List[str] = field(default_factory=list)
    context_notes: List[str] = field(default_factory=list)
    biological_context: Dict[str, Any] = field(default_factory=dict)

    # Phase 4A Explicit Five-Dimension Evidence Segregation
    @property
    def provenance_validity(self) -> str:
        return self.status.value

    @property
    def context_compatibility(self) -> str:
        return self.context_status.value

    @property
    def data_availability(self) -> str:
        if self.status in [ValidationStatusCode.FILE_NOT_FOUND, ValidationStatusCode.EXPERIMENT_NOT_FOUND, ValidationStatusCode.FILE_UNAVAILABLE, ValidationStatusCode.RATE_LIMIT_EXCEEDED, ValidationStatusCode.SERVER_ERROR]:
            return "DATA_UNAVAILABLE"
        return "DATA_AVAILABLE"

    def to_evidence_summary(self, peak_detected: Optional[bool] = None) -> Dict[str, str]:
        """Returns the 5 separate unaggregated scientific evaluation dimensions."""
        obs = "PEAK_DETECTED" if peak_detected is True else ("ASSAYED_NO_PEAK_DETECTED" if peak_detected is False else "NOT_EVALUATED")
        interp = "PROXIMAL_OCCUPANCY_OBSERVED" if peak_detected is True else ("NO_PROXIMAL_OCCUPANCY_IN_EVALUATED_WINDOW" if peak_detected is False else "INSUFFICIENT_DATA")
        return {
            "PROVENANCE_VALIDITY": self.provenance_validity,
            "CONTEXT_COMPATIBILITY": self.context_compatibility,
            "DATA_AVAILABILITY": self.data_availability,
            "EXPERIMENTAL_OBSERVATION": obs,
            "REGULATORY_INTERPRETATION": interp
        }

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["status"] = self.status.value
        d["checksum_status"] = self.checksum_status.value
        d["context_status"] = self.context_status.value
        d["provenance_validity"] = self.provenance_validity
        d["context_compatibility"] = self.context_compatibility
        d["data_availability"] = self.data_availability
        return d



class EvidenceIntegrityEngine:
    """
    Automated scientific evidence verification engine.
    Cross-checks user-supplied ChIP-seq datasets, files, and experiments
    against authoritative ENCODE metadata schemas with multi-tier caching and rate-limit handling.
    """

    def __init__(self, cache_dir: Optional[str] = None, memory_snapshot: Optional[Dict[str, Any]] = None):
        self.cache_dir = cache_dir or os.path.join(os.getcwd(), ".tf_explorer_cache", "encode_metadata")
        os.makedirs(self.cache_dir, exist_ok=True)
        self.memory_snapshot = memory_snapshot if memory_snapshot is not None else VERIFIED_METADATA_SNAPSHOT

    def _normalize_tf(self, tf: str) -> str:
        t = tf.strip().upper()
        return TF_ALIASES.get(t, t)

    def _normalize_assembly(self, assembly: str) -> str:
        a = assembly.strip().upper()
        return ASSEMBLY_ALIASES.get(a, assembly.strip())

    def _normalize_biosample(self, biosample: str) -> str:
        return biosample.strip().replace(" ", "").replace("-", "").upper()

    def _safe_encode_get(self, url: str, max_retries: int = 3) -> Tuple[Optional[requests.Response], Optional[str]]:
        """
        Executes HTTP GET against ENCODE with exponential backoff on HTTP 429 and 5xx responses.
        Returns (response, error_code).
        """
        headers = {"Accept": "application/json", "User-Agent": "TF-Explorer/3.2 (Academic; Integrity Engine)"}
        for attempt in range(max_retries):
            try:
                resp = requests.get(url, headers=headers, timeout=10)
                if resp.status_code == 200:
                    return resp, None
                elif resp.status_code == 404:
                    return resp, "404_NOT_FOUND"
                elif resp.status_code == 429:
                    retry_after = float(resp.headers.get("Retry-After", 1.5 * (attempt + 1)))
                    logger.warning(f"ENCODE rate limit encountered on {url}. Retrying in {retry_after}s...")
                    time.sleep(retry_after)
                    if attempt == max_retries - 1:
                        return resp, "429_RATE_LIMIT"
                elif resp.status_code >= 500:
                    logger.warning(f"ENCODE server error {resp.status_code} on {url}. Attempt {attempt+1}/{max_retries}...")
                    time.sleep(1.0 * (attempt + 1))
                    if attempt == max_retries - 1:
                        return resp, f"SERVER_ERROR_{resp.status_code}"
                else:
                    return resp, f"HTTP_{resp.status_code}"
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
                logger.warning(f"Connection error to ENCODE: {e}. Attempt {attempt+1}/{max_retries}...")
                time.sleep(1.0 * (attempt + 1))
                if attempt == max_retries - 1:
                    return None, f"CONNECTION_ERROR: {str(e)}"
            except Exception as e:
                return None, f"REQUEST_ERROR: {str(e)}"
        return None, "MAX_RETRIES_EXCEEDED"

    def get_file_metadata(self, file_accession: str, offline_only: bool = False) -> Tuple[Optional[Dict[str, Any]], str]:
        """
        Retrieves authoritative ENCODE file metadata.
        Checks memory snapshot, then local disk cache, then live ENCODE REST API.
        """
        acc = file_accession.strip()
        v, msg = validate_accession_format(acc, "encode_file")
        if not v:
            return None, f"Invalid file accession format: {acc}"

        # 1. Check memory snapshot
        known_files = self.memory_snapshot.get("encode_files", {})
        if acc in known_files:
            rec = dict(known_files[acc])
            rec.setdefault("file_accession", acc)
            return rec, "memory_snapshot"

        # 2. Check disk cache
        cache_file = os.path.join(self.cache_dir, f"file_{acc}.json")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached_data = json.load(f)
                    # Verify cache schema integrity
                    if isinstance(cached_data, dict) and cached_data.get("file_accession") == acc:
                        return cached_data, "disk_cache"
            except Exception as e:
                logger.warning(f"Error reading cached file metadata for {acc}: {e}")

        if offline_only:
            return None, "File not present in offline cache/snapshot"

        # 3. Live query from ENCODE REST API
        url = f"{ENCODE_BASE_URL}/files/{acc}/?format=json"
        resp, err_code = self._safe_encode_get(url)

        if resp is None:
            return None, f"Live query failed: {err_code}"
        if resp.status_code == 404:
            return None, f"File {acc} not found in ENCODE"
        if resp.status_code != 200:
            return None, f"ENCODE query failed with status {resp.status_code}: {err_code}"

        try:
            raw_text = resp.text
            raw_data = resp.json()

            dataset_field = raw_data.get("dataset", "")
            exp_acc = dataset_field.strip("/").split("/")[-1] if dataset_field else ""

            target_label = raw_data.get("target", {}).get("label") if isinstance(raw_data.get("target"), dict) else None
            biosample_term = raw_data.get("biosample_ontology", {}).get("term_name") if isinstance(raw_data.get("biosample_ontology"), dict) else None

            record = {
                "cache_schema_version": CACHE_SCHEMA_VERSION,
                "file_accession": acc,
                "source_url": url,
                "retrieval_timestamp": datetime.now(timezone.utc).isoformat(),
                "content_hash": hashlib.sha256(raw_text.encode("utf-8")).hexdigest(),
                "experiment": exp_acc,
                "target": target_label,
                "biosample": biosample_term,
                "assembly": raw_data.get("assembly", ""),
                "format": raw_data.get("file_format", ""),
                "format_type": raw_data.get("file_format_type", ""),
                "output_type": raw_data.get("output_type", ""),
                "status": raw_data.get("status", ""),
                "md5sum": raw_data.get("md5sum", ""),
                "assay_title": raw_data.get("assay_title", ""),
                "raw_dataset": dataset_field
            }

            # Cache to disk
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(record, f, indent=2)

            return record, "encode_rest_api"
        except Exception as e:
            logger.warning(f"Error parsing live response for file {acc}: {e}")
            return None, f"Parsing error: {str(e)}"

    def get_experiment_metadata(self, experiment_accession: str, offline_only: bool = False) -> Tuple[Optional[Dict[str, Any]], str]:
        """
        Retrieves authoritative ENCODE experiment metadata.
        """
        acc = experiment_accession.strip()
        v, msg = validate_accession_format(acc, "encode_experiment")
        if not v:
            return None, f"Invalid experiment accession format: {acc}"

        # 1. Check memory snapshot
        known_exps = self.memory_snapshot.get("encode_experiments", {})
        if acc in known_exps:
            raw_rec = known_exps[acc]
            treatment_metadata_present = ("is_treated" in raw_rec) or ("treatments" in raw_rec)
            gm_metadata_present = ("is_genetically_modified" in raw_rec) or ("genetic_modifications" in raw_rec)
            rec = dict(raw_rec)
            rec.setdefault("experiment_accession", acc)
            raw_tr = rec.get("treatments", [])
            tr_names = []
            tr_detailed = []
            for t in raw_tr:
                if isinstance(t, dict):
                    tr_detailed.append(t)
                    t_name = t.get("treatment_term_name", "")
                    if t_name and t_name not in tr_names:
                        tr_names.append(t_name)
                elif isinstance(t, str):
                    if t and t not in tr_names:
                        tr_names.append(t)
                    tr_detailed.append({"treatment_term_name": t})
            rec["treatments"] = tr_names
            rec["treatments_detailed"] = rec.get("treatments_detailed") or tr_detailed
            rec["is_treated"] = len(tr_names) > 0 or rec.get("is_treated", False)
            rec.setdefault("is_genetically_modified", len(rec.get("genetic_modifications", [])) > 0)
            rec.setdefault("genetic_modifications", [])
            rec.setdefault("biosample_summary", "")
            rec.setdefault("biosample_ontology_id", "")
            rec["treatment_metadata_present"] = rec.get("treatment_metadata_present", treatment_metadata_present)
            rec["gm_metadata_present"] = rec.get("gm_metadata_present", gm_metadata_present)
            return rec, "memory_snapshot"

        # 2. Check disk cache
        cache_file = os.path.join(self.cache_dir, f"exp_{acc}.json")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached_data = json.load(f)
                    if isinstance(cached_data, dict) and cached_data.get("experiment_accession") == acc:
                        treatment_metadata_present = ("is_treated" in cached_data) or ("treatments" in cached_data)
                        gm_metadata_present = ("is_genetically_modified" in cached_data) or ("genetic_modifications" in cached_data)
                        raw_tr = cached_data.get("treatments", [])
                        tr_names = []
                        tr_detailed = []
                        for t in raw_tr:
                            if isinstance(t, dict):
                                tr_detailed.append(t)
                                t_name = t.get("treatment_term_name", "")
                                if t_name and t_name not in tr_names:
                                    tr_names.append(t_name)
                            elif isinstance(t, str):
                                if t and t not in tr_names:
                                    tr_names.append(t)
                                tr_detailed.append({"treatment_term_name": t})
                        cached_data["treatments"] = tr_names
                        cached_data["treatments_detailed"] = cached_data.get("treatments_detailed") or tr_detailed
                        cached_data["is_treated"] = len(tr_names) > 0 or cached_data.get("is_treated", False)
                        cached_data.setdefault("is_genetically_modified", len(cached_data.get("genetic_modifications", [])) > 0)
                        cached_data.setdefault("genetic_modifications", [])
                        cached_data.setdefault("replicates", [])
                        cached_data.setdefault("biosample_summary", "")
                        cached_data.setdefault("biosample_ontology_id", "")
                        cached_data["treatment_metadata_present"] = cached_data.get("treatment_metadata_present", treatment_metadata_present)
                        cached_data["gm_metadata_present"] = cached_data.get("gm_metadata_present", gm_metadata_present)
                        return cached_data, "disk_cache"
            except Exception as e:
                logger.warning(f"Error reading cached experiment metadata for {acc}: {e}")



        if offline_only:
            return None, "Experiment not present in offline cache/snapshot"

        # 3. Live query from ENCODE REST API
        url = f"{ENCODE_BASE_URL}/experiments/{acc}/?format=json"
        resp, err_code = self._safe_encode_get(url)

        if resp is None:
            return None, f"Live query failed: {err_code}"
        if resp.status_code == 404:
            return None, f"Experiment {acc} not found in ENCODE"
        if resp.status_code != 200:
            return None, f"ENCODE query failed with status {resp.status_code}: {err_code}"

        try:
            raw_text = resp.text
            raw_data = resp.json()

            target_label = raw_data.get("target", {}).get("label") if isinstance(raw_data.get("target"), dict) else "Unknown"
            biosample_term = raw_data.get("biosample_ontology", {}).get("term_name") if isinstance(raw_data.get("biosample_ontology"), dict) else "Unknown"
            biosample_onto_id = raw_data.get("biosample_ontology", {}).get("term_id", "") if isinstance(raw_data.get("biosample_ontology"), dict) else ""
            biosample_summary = raw_data.get("biosample_summary", "")

            # Extract detailed treatments and genetic modifications across replicates
            treatments_names = []
            treatments_detailed = []
            genetic_modifications_detailed = []

            replicates_summary = []
            for rep in raw_data.get("replicates", []):
                rep_num = rep.get("biological_replicate_number", 1)
                lib = rep.get("library", {})
                bio = lib.get("biosample", {})
                rep_tr = []
                for tr in bio.get("treatments", []):
                    t_name = tr.get("treatment_term_name", "")
                    if t_name:
                        rep_tr.append(t_name)
                    if t_name and t_name not in treatments_names:
                        treatments_names.append(t_name)
                    t_rec = {
                        "treatment_term_name": t_name,
                        "treatment_term_id": tr.get("treatment_term_id", ""),
                        "amount": tr.get("amount"),
                        "amount_units": tr.get("amount_units", ""),
                        "duration": tr.get("duration"),
                        "duration_units": tr.get("duration_units", "")
                    }
                    if t_rec not in treatments_detailed:
                        treatments_detailed.append(t_rec)

                rep_gms = []
                for gm in bio.get("genetic_modifications", []):
                    tags = [tag.get("name", "tag") for tag in gm.get("introduced_tags", [])]
                    if tags:
                        rep_gms.extend(tags)
                    elif gm.get("category"):
                        rep_gms.append(gm.get("category"))
                    gm_rec = {
                        "category": gm.get("category", ""),
                        "purpose": gm.get("purpose", ""),
                        "method": gm.get("method", ""),
                        "tags": gm.get("introduced_tags", [])
                    }
                    if gm_rec not in genetic_modifications_detailed:
                        genetic_modifications_detailed.append(gm_rec)

                replicates_summary.append({
                    "biological_replicate_number": rep_num,
                    "treatments": rep_tr,
                    "genetic_modifications": rep_gms,
                    "is_treated": len(rep_tr) > 0,
                    "is_genetically_modified": len(rep_gms) > 0
                })

            is_treated = len(treatments_detailed) > 0
            is_gm = len(genetic_modifications_detailed) > 0

            files_list = []
            for f in raw_data.get("files", []):
                if isinstance(f, dict):
                    files_list.append(f.get("accession"))
                elif isinstance(f, str):
                    files_list.append(f.strip("/").split("/")[-1])

            record = {
                "cache_schema_version": CACHE_SCHEMA_VERSION,
                "experiment_accession": acc,
                "source_url": url,
                "retrieval_timestamp": datetime.now(timezone.utc).isoformat(),
                "content_hash": hashlib.sha256(raw_text.encode("utf-8")).hexdigest(),
                "target": target_label,
                "biosample": biosample_term,
                "biosample_ontology_id": biosample_onto_id,
                "biosample_summary": biosample_summary,
                "treatments": treatments_names,
                "treatments_detailed": treatments_detailed,
                "genetic_modifications": genetic_modifications_detailed,
                "replicates": replicates_summary,
                "is_treated": is_treated,
                "is_genetically_modified": is_gm,
                "assay_title": raw_data.get("assay_title", "TF ChIP-seq"),
                "status": raw_data.get("status", "released"),
                "description": raw_data.get("description", ""),
                "files": files_list
            }

            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(record, f, indent=2)

            return record, "encode_rest_api"
        except Exception as e:
            logger.warning(f"Error parsing live response for experiment {acc}: {e}")
            return None, f"Parsing error: {str(e)}"

    def verify(
        self,
        file_accession: str,
        experiment_accession: str,
        expected_tf: str,
        expected_biosample: str,
        expected_assembly: str = "GRCh38",
        expected_assay: Optional[str] = "TF ChIP-seq",
        expected_output_type: Optional[str] = None,
        local_file_path: Optional[str] = None,
        expected_md5sum: Optional[str] = None,
        offline_only: bool = False,
        treatment_requirement: Optional[str] = None,
        genetic_modification_requirement: Optional[str] = None,
        allow_exploratory_context: bool = False,
        expected_untreated: Optional[bool] = None,
        expected_unmodified: Optional[bool] = None
    ) -> ValidationResult:
        """
        Executes exhaustive 12-point integrity check on dataset provenance.
        """
        if expected_untreated is True and treatment_requirement is None:
            treatment_requirement = "untreated"
        if expected_unmodified is True and genetic_modification_requirement is None:
            genetic_modification_requirement = "unmodified"

        accessions_dict = {"file": file_accession, "experiment": experiment_accession}

        # Check 1: Accession format
        vf, msg_f = validate_accession_format(file_accession, "encode_file")
        if not vf:
            return ValidationResult(
                status=ValidationStatusCode.INVALID_ACCESSION_FORMAT,
                is_valid=False,
                failure_code="INVALID_ACCESSION_FORMAT",
                checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                accessions=accessions_dict,
                message=f"File accession syntax rejected: {msg_f}"
            )

        ve, msg_e = validate_accession_format(experiment_accession, "encode_experiment")
        if not ve:
            return ValidationResult(
                status=ValidationStatusCode.INVALID_ACCESSION_FORMAT,
                is_valid=False,
                failure_code="INVALID_ACCESSION_FORMAT",
                checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                accessions=accessions_dict,
                message=f"Experiment accession syntax rejected: {msg_e}"
            )

        # Check 2: Fetch File Metadata
        f_meta, f_source = self.get_file_metadata(file_accession, offline_only=offline_only)
        if not f_meta:
            fail_status = ValidationStatusCode.FILE_NOT_FOUND
            if "RATE_LIMIT" in f_source:
                fail_status = ValidationStatusCode.RATE_LIMIT_EXCEEDED
            elif "SERVER_ERROR" in f_source:
                fail_status = ValidationStatusCode.SERVER_ERROR
            elif offline_only:
                fail_status = ValidationStatusCode.METADATA_UNVERIFIED

            return ValidationResult(
                status=fail_status,
                is_valid=False,
                failure_code=fail_status.value,
                checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                accessions=accessions_dict,
                message=f"File {file_accession} could not be resolved from authoritative metadata: {f_source}"
            )

        # Check 3: Fetch Experiment Metadata
        e_meta, e_source = self.get_experiment_metadata(experiment_accession, offline_only=offline_only)
        if not e_meta:
            fail_status = ValidationStatusCode.EXPERIMENT_NOT_FOUND
            if "RATE_LIMIT" in e_source:
                fail_status = ValidationStatusCode.RATE_LIMIT_EXCEEDED
            elif "SERVER_ERROR" in e_source:
                fail_status = ValidationStatusCode.SERVER_ERROR
            elif offline_only:
                fail_status = ValidationStatusCode.METADATA_UNVERIFIED

            return ValidationResult(
                status=fail_status,
                is_valid=False,
                failure_code=fail_status.value,
                checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                accessions=accessions_dict,
                message=f"Experiment {experiment_accession} could not be resolved: {e_source}"
            )

        meta_source = f"{f_source} + {e_source}"
        evidence_dict = {
            "file": f_meta,
            "experiment": e_meta
        }

        # Check 4: File belongs to Parent Experiment
        f_parent = f_meta.get("experiment", "")
        if f_parent and f_parent != experiment_accession:
            return ValidationResult(
                status=ValidationStatusCode.FILE_EXPERIMENT_MISMATCH,
                is_valid=False,
                failure_code="FILE_EXPERIMENT_MISMATCH",
                checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                metadata_source=meta_source,
                accessions=accessions_dict,
                message=f"File '{file_accession}' belongs to parent experiment '{f_parent}', not '{experiment_accession}'.",
                evidence=evidence_dict
            )

        # Check 5: Target TF Matching
        norm_expected_tf = self._normalize_tf(expected_tf)
        exp_target = e_meta.get("target") or f_meta.get("target") or ""
        norm_actual_tf = self._normalize_tf(exp_target)

        if norm_actual_tf != norm_expected_tf:
            return ValidationResult(
                status=ValidationStatusCode.TF_TARGET_MISMATCH,
                is_valid=False,
                failure_code="TF_TARGET_MISMATCH",
                checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                metadata_source=meta_source,
                accessions=accessions_dict,
                message=f"Experiment '{experiment_accession}' targets TF '{exp_target}', which does not match requested '{expected_tf}'.",
                evidence=evidence_dict
            )

        # Check 6: Biosample Matching
        norm_expected_bio = self._normalize_biosample(expected_biosample)
        actual_bio = e_meta.get("biosample") or f_meta.get("biosample") or ""
        norm_actual_bio = self._normalize_biosample(actual_bio)

        if norm_expected_bio not in norm_actual_bio and norm_actual_bio not in norm_expected_bio:
            return ValidationResult(
                status=ValidationStatusCode.BIOSAMPLE_MISMATCH,
                is_valid=False,
                failure_code="BIOSAMPLE_MISMATCH",
                checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                metadata_source=meta_source,
                accessions=accessions_dict,
                message=f"Experiment '{experiment_accession}' was assayed in biosample '{actual_bio}', not '{expected_biosample}'.",
                evidence=evidence_dict
            )

        # Check 7: Assembly Compatibility
        norm_expected_asm = self._normalize_assembly(expected_assembly)
        actual_asm = f_meta.get("assembly", "")
        norm_actual_asm = self._normalize_assembly(actual_asm) if actual_asm else norm_expected_asm

        if norm_actual_asm != norm_expected_asm:
            return ValidationResult(
                status=ValidationStatusCode.ASSEMBLY_MISMATCH,
                is_valid=False,
                failure_code="ASSEMBLY_MISMATCH",
                checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                metadata_source=meta_source,
                accessions=accessions_dict,
                message=f"File '{file_accession}' is mapped to genome assembly '{actual_asm}', cannot be evaluated on '{expected_assembly}'.",
                evidence=evidence_dict
            )

        # Check 8: Assay Type Matching
        if expected_assay:
            actual_assay = e_meta.get("assay_title") or f_meta.get("assay_title") or ""
            if actual_assay and actual_assay.upper() != expected_assay.upper():
                return ValidationResult(
                    status=ValidationStatusCode.ASSAY_TYPE_MISMATCH,
                    is_valid=False,
                    failure_code="ASSAY_TYPE_MISMATCH",
                    checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                    metadata_source=meta_source,
                    accessions=accessions_dict,
                    message=f"Experiment assay type '{actual_assay}' does not match expected '{expected_assay}'.",
                    evidence=evidence_dict
                )

        # Check 9: Output Type Matching
        if expected_output_type:
            actual_output = f_meta.get("output_type", "")
            if actual_output and actual_output.lower() != expected_output_type.lower():
                return ValidationResult(
                    status=ValidationStatusCode.OUTPUT_TYPE_MISMATCH,
                    is_valid=False,
                    failure_code="OUTPUT_TYPE_MISMATCH",
                    checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                    metadata_source=meta_source,
                    accessions=accessions_dict,
                    message=f"File output type '{actual_output}' does not match expected '{expected_output_type}'.",
                    evidence=evidence_dict
                )

        # Check 10: Status check
        f_status = f_meta.get("status", "released").lower()
        if f_status in ("revoked", "replaced", "deleted"):
            return ValidationResult(
                status=ValidationStatusCode.STATUS_DEPRECATED_OR_REVOKED,
                is_valid=False,
                failure_code="STATUS_DEPRECATED_OR_REVOKED",
                checksum_status=ChecksumStatusCode.FILE_NOT_DOWNLOADED,
                metadata_source=meta_source,
                accessions=accessions_dict,
                message=f"File '{file_accession}' has deprecated/invalid processing status: '{f_status}'.",
                evidence=evidence_dict
            )

        # Check 11: MD5 Checksum Verification (Explicit Separation)
        auth_md5 = f_meta.get("md5sum", "").lower()
        chk_status = ChecksumStatusCode.FILE_NOT_DOWNLOADED

        if local_file_path and os.path.exists(local_file_path):
            hasher = hashlib.md5()
            with open(local_file_path, "rb") as lf:
                for chunk in iter(lambda: lf.read(65536), b""):
                    hasher.update(chunk)
            local_md5 = hasher.hexdigest().lower()
            if auth_md5 and local_md5 != auth_md5:
                return ValidationResult(
                    status=ValidationStatusCode.CHECKSUM_MISMATCH,
                    is_valid=False,
                    failure_code="CHECKSUM_MISMATCH",
                    checksum_status=ChecksumStatusCode.CHECKSUM_MISMATCH,
                    metadata_source=meta_source,
                    accessions=accessions_dict,
                    message=f"Local file '{local_file_path}' MD5 '{local_md5}' does not match ENCODE checksum '{auth_md5}'.",
                    evidence=evidence_dict
                )
            chk_status = ChecksumStatusCode.LOCAL_FILE_CHECKSUM_VERIFIED
        elif expected_md5sum:
            if auth_md5 and expected_md5sum.lower() != auth_md5:
                return ValidationResult(
                    status=ValidationStatusCode.CHECKSUM_MISMATCH,
                    is_valid=False,
                    failure_code="CHECKSUM_MISMATCH",
                    checksum_status=ChecksumStatusCode.CHECKSUM_MISMATCH,
                    metadata_source=meta_source,
                    accessions=accessions_dict,
                    message=f"Expected MD5 '{expected_md5sum}' does not match authoritative checksum '{auth_md5}'.",
                    evidence=evidence_dict
                )
            chk_status = ChecksumStatusCode.REMOTE_CHECKSUM_RETRIEVED
        else:
            chk_status = ChecksumStatusCode.FILE_NOT_DOWNLOADED

        # Check 12: Experimental Context Evaluation (Phase 3C-A / Phase 3C-B)
        actual_is_treated = e_meta.get("is_treated", len(e_meta.get("treatments", [])) > 0)
        actual_is_gm = e_meta.get("is_genetically_modified", len(e_meta.get("genetic_modifications", [])) > 0)
        actual_treatments = e_meta.get("treatments_detailed", e_meta.get("treatments", []))
        actual_gms = e_meta.get("genetic_modifications", [])
        actual_bio_summary = e_meta.get("biosample_summary", "")
        actual_bio_onto_id = e_meta.get("biosample_ontology_id", "")

        context_warnings = []
        context_mismatch_reasons = []

        # Check metadata completeness
        has_treatment_meta = e_meta.get("treatment_metadata_present", ("is_treated" in e_meta) or ("treatments" in e_meta))
        has_gm_meta = e_meta.get("gm_metadata_present", ("is_genetically_modified" in e_meta) or ("genetic_modifications" in e_meta))
        metadata_incomplete_reasons = []

        if treatment_requirement is not None and not has_treatment_meta:
            metadata_incomplete_reasons.append("Experiment lacks authoritative treatment metadata records.")
        if genetic_modification_requirement is not None and not has_gm_meta:
            metadata_incomplete_reasons.append("Experiment lacks authoritative genetic modification metadata records.")

        # Check replicate discordance
        reps = e_meta.get("replicates", [])
        has_replicate_discordance = False
        if len(reps) > 1:
            tr_states = [r.get("is_treated") for r in reps if "is_treated" in r]
            gm_states = [r.get("is_genetically_modified") for r in reps if "is_genetically_modified" in r]
            if len(set(tr_states)) > 1 or len(set(gm_states)) > 1:
                has_replicate_discordance = True
                context_warnings.append(
                    "⚠️ REPLICATE CONTEXT DISCORDANCE: Biological replicates have discordant treatment or modification annotations. "
                    "Findings must not assume uniform experimental conditions across replicates."
                )

        if treatment_requirement is not None and has_treatment_meta:
            norm_tr = treatment_requirement.strip().upper()
            if norm_tr in ("UNTREATED", "NONE", "FALSE"):
                if actual_is_treated:
                    tr_names = [t.get("treatment_term_name", str(t)) if isinstance(t, dict) else str(t) for t in actual_treatments]
                    context_mismatch_reasons.append(
                        f"Sample was treated with: {', '.join(tr_names)} (expected untreated)."
                    )
            elif norm_tr not in ("ANY", "ALL"):
                tr_names = [t.get("treatment_term_name", str(t)).upper() if isinstance(t, dict) else str(t).upper() for t in actual_treatments]
                if not any(norm_tr in tn for tn in tr_names):
                    context_mismatch_reasons.append(
                        f"Requested treatment '{treatment_requirement}' not found in actual treatments: {tr_names}."
                    )

        if genetic_modification_requirement is not None and has_gm_meta:
            norm_gm = genetic_modification_requirement.strip().upper()
            if norm_gm in ("UNMODIFIED", "WILDTYPE", "WILD_TYPE", "WT", "FALSE"):
                if actual_is_gm:
                    gm_tags = []
                    for gm in actual_gms:
                        if isinstance(gm, dict):
                            for tag in gm.get("tags", []):
                                gm_tags.append(tag.get("name", "tag"))
                    tag_str = f" (tags: {', '.join(gm_tags)})" if gm_tags else ""
                    context_mismatch_reasons.append(
                        f"Sample is genetically modified{tag_str} (expected unmodified)."
                    )
            elif norm_gm not in ("ANY", "ALL"):
                gm_tags = []
                for gm in actual_gms:
                    if isinstance(gm, dict):
                        for tag in gm.get("tags", []):
                            gm_tags.append(tag.get("name", "").upper())
                if not any(norm_gm in gt for gt in gm_tags):
                    context_mismatch_reasons.append(
                        f"Requested tag/modification '{genetic_modification_requirement}' not found in sample: {gm_tags}."
                    )

        if metadata_incomplete_reasons:
            if allow_exploratory_context:
                context_status = ExperimentalContextStatus.CONTEXT_ACCEPTED_WITH_LIMITATIONS
                context_compatible = True
                context_accepted = True
                context_warnings.extend(metadata_incomplete_reasons)
                context_warnings.append(
                    "EXPLORATORY CONTEXT QUALIFICATION: Analysis accepted despite incomplete experimental context metadata."
                )
            else:
                context_status = ExperimentalContextStatus.CONTEXT_METADATA_INCOMPLETE
                context_compatible = False
                context_accepted = False
                context_warnings.extend(metadata_incomplete_reasons)
        elif treatment_requirement is None and genetic_modification_requirement is None:
            context_status = ExperimentalContextStatus.CONTEXT_UNSPECIFIED
            context_compatible = True
            context_accepted = False
            if actual_is_treated:
                context_warnings.append(f"Note: Experiment utilizes treated biosample ({actual_bio_summary}).")
            if actual_is_gm:
                context_warnings.append(f"Note: Experiment utilizes genetically modified/tagged biosample ({actual_bio_summary}).")
        elif context_mismatch_reasons:
            if allow_exploratory_context:
                context_status = ExperimentalContextStatus.CONTEXT_ACCEPTED_WITH_LIMITATIONS
                context_compatible = True
                context_accepted = True
                context_warnings.extend(context_mismatch_reasons)
                context_warnings.append(
                    "EXPLORATORY CONTEXT QUALIFICATION: Analysis accepted despite experimental context divergence. "
                    "Findings must not be attributed to basal untreated/unmodified conditions."
                )
            else:
                context_status = ExperimentalContextStatus.CONTEXT_MISMATCH
                context_compatible = False
                context_accepted = False
                context_warnings.extend(context_mismatch_reasons)
        else:
            context_status = ExperimentalContextStatus.CONTEXT_MATCH
            context_compatible = True
            context_accepted = False

        biological_context = {
            "biosample_term": actual_bio,
            "biosample_ontology_id": actual_bio_onto_id,
            "biosample_summary": actual_bio_summary,
            "is_treated": actual_is_treated,
            "treatments": actual_treatments,
            "is_genetically_modified": actual_is_gm,
            "genetic_modifications": actual_gms,
            "has_replicate_discordance": has_replicate_discordance,
            "replicates": reps
        }

        # Check 13: All criteria satisfied
        return ValidationResult(
            status=ValidationStatusCode.VALIDATED,
            is_valid=True,
            checksum_status=chk_status,
            metadata_source=meta_source,
            accessions=accessions_dict,
            message=f"Dataset verified authentic: {file_accession} ({norm_expected_tf} in {actual_bio}, {norm_expected_asm}).",
            evidence=evidence_dict,
            context_status=context_status,
            context_compatible=context_compatible,
            context_accepted=context_accepted,
            context_warnings=context_warnings,
            context_notes=context_warnings,
            biological_context=biological_context
        )

