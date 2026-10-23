"""
Validated Configuration Factory for TF-Explorer Phase 4B-R.

Constructs fully authenticated, source-backed EvidenceConfiguration objects.
Guarantees that observation status, peak overlap, promoter coordinates, and
evidence interpretations are derived directly from authoritative primary data
rather than preassigned or guessed by callers.
"""

import os
import gzip
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Union

from tf_explorer.disagreement_taxonomy import (
    EvidenceConfiguration,
    EvidenceInterpretation,
    calculate_signed_promoter_coordinates
)
from tf_explorer.disagreement_provenance import (
    verify_source_identity_contract,
    create_synthetic_provenance_token_for_testing,
    scan_bed_for_promoter_overlap,
    resolve_transcript_version,
    AuthoritativeIdentityStatus,
    BiologicalRelationshipStatus,
    LocalParsingStatus,
    ProvenanceToken,
    AUTHORITATIVE_ENCODE_CATALOG
)
from tf_explorer.peak_authenticator import PrimaryPeakAuthenticator

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def find_source_bed_file(
    file_accession: str,
    search_dirs: Optional[List[Path]] = None
) -> Optional[Path]:
    """Locates the local primary BED file for a given ENCODE file accession."""
    if search_dirs is None:
        search_dirs = [
            PROJECT_ROOT / "encode_cache",
            PROJECT_ROOT / "reference_data" / "peaks",
            PROJECT_ROOT / "reference_data",
        ]

    for d in search_dirs:
        if not d.exists():
            continue
        candidates = [
            d / f"{file_accession}.bed.gz",
            d / f"{file_accession}.bed",
            d / f"{file_accession}.narrowPeak.gz",
            d / f"{file_accession}.narrowPeak",
            d / file_accession
        ]
        for c in candidates:
            if c.is_file() and c.stat().st_size > 0:
                return c
    return None


def verify_source_bed_integrity(
    bed_file: Path,
    authenticator: PrimaryPeakAuthenticator
) -> Tuple[bool, Dict[str, Any]]:
    """Verifies file readability, cryptographic hashes, and gzip stream integrity."""
    if not bed_file.exists() or bed_file.stat().st_size == 0:
        return False, {"message": "File does not exist or is empty"}

    try:
        md5_val, sha256_val = authenticator.compute_hashes(bed_file)
        file_acc = bed_file.name.split(".")[0]
        # Validate that gzip/BED can be parsed without error
        authenticator.load_peaks(file_acc, bed_file)
        return True, {"md5": md5_val, "sha256": sha256_val}
    except Exception as e:
        logger.warning(f"Integrity check failed for {bed_file}: {e}")
        return False, {"message": str(e)}


def build_validated_configuration(
    gene_symbol: str,
    transcript_id: str,
    canonical_tss: int,
    strand: str,
    tf_target: str,
    encode_experiment_accession: str,
    encode_file_accession: str,
    biosample: str,
    treatment_context: Optional[str] = None,
    genetic_modification_context: Optional[str] = None,
    processing_tier: str = "conservative_idr_peaks",
    window_name: str = "Standard",
    window_start_offset: int = -1000,
    window_end_offset: int = 200,
    chromosome: str = "",
    config_id: Optional[str] = None,
    authenticator: Optional[PrimaryPeakAuthenticator] = None,
    assembly: str = "GRCh38",
    is_mane_select: bool = False,
    is_mane_plus_clinical: bool = False,
    mane_status: str = "MANE_UNVERIFIED",
    custom_bed_path: Optional[Path] = None
) -> EvidenceConfiguration:
    """
    Factory constructing a fully verified, source-backed EvidenceConfiguration.
    Calculates promoter bounds, authenticates source files, scans primary peak records,
    and derives observation status and evidence interpretations fail-closed.
    """
    meta = AUTHORITATIVE_ENCODE_CATALOG.get(encode_file_accession.strip().upper())
    explicit_treatment = treatment_context is not None
    explicit_modification = genetic_modification_context is not None

    if treatment_context is None:
        treatment_context = meta.get("treatment_context", "Untreated") if meta else "Untreated"
    if genetic_modification_context is None:
        genetic_modification_context = meta.get("genetic_modification_context", "None") if meta else "None"

    if assembly != "GRCh38":
        p_s, p_e = calculate_signed_promoter_coordinates(
            canonical_tss,
            strand if strand in ["+", "-"] else "+",
            window_start_offset,
            window_end_offset,
            is_bed0=True
        )
        return EvidenceConfiguration(
            config_id=config_id or f"{gene_symbol}_{biosample}_{window_name}_{encode_file_accession}",
            gene_symbol=gene_symbol.upper(),
            transcript_id=transcript_id,
            canonical_tss=canonical_tss,
            strand=strand,
            assembly=assembly,
            chromosome=chromosome,
            is_mane_select=is_mane_select,
            is_mane_plus_clinical=is_mane_plus_clinical,
            mane_status=mane_status,
            window_name=window_name,
            window_start_offset=window_start_offset,
            window_end_offset=window_end_offset,
            promoter_upstream=abs(window_start_offset),
            promoter_downstream=abs(window_end_offset),
            promoter_start=p_s,
            promoter_end=p_e,
            tf_target=tf_target.upper(),
            biosample=biosample,
            encode_experiment_accession=encode_experiment_accession,
            encode_file_accession=encode_file_accession,
            treatment_context=treatment_context,
            genetic_modification_context=genetic_modification_context,
            processing_tier=processing_tier,
            is_authenticated=False,
            auth_status="ASSEMBLY_MISMATCH",
            observation_status="SOURCE_NOT_VERIFIED",
            evidence_status=EvidenceInterpretation.SOURCE_NOT_VERIFIED,
            exclusion_reason=f"Unsupported genome assembly '{assembly}'. TF-Explorer strictly requires GRCh38."
        )

    if strand not in ["+", "-"]:
        raise ValueError(f"Strand must be '+' or '-'. Got: {strand}")

    if window_start_offset >= window_end_offset:
        raise ValueError(
            f"window_start_offset ({window_start_offset}) must be strictly less than "
            f"window_end_offset ({window_end_offset})"
        )

    # 1. Calculate BED0 promoter genomic bounds
    promoter_start, promoter_end = calculate_signed_promoter_coordinates(
        canonical_tss,
        strand,
        window_start_offset,
        window_end_offset,
        is_bed0=True
    )

    if config_id is None:
        clean_win = window_name.replace(" ", "_").replace("[", "").replace("]", "")
        config_id = f"{gene_symbol}_{biosample}_{clean_win}_{encode_file_accession}"

    # 2. Locate source BED file
    bed_file = custom_bed_path or find_source_bed_file(encode_file_accession)
    if authenticator is None:
        authenticator = PrimaryPeakAuthenticator()

    # Fail closed if file is missing
    if bed_file is None or not bed_file.exists():
        return EvidenceConfiguration(
            config_id=config_id,
            gene_symbol=gene_symbol.upper(),
            transcript_id=transcript_id,
            canonical_tss=canonical_tss,
            strand=strand,
            assembly=assembly,
            chromosome=chromosome,
            is_mane_select=is_mane_select,
            is_mane_plus_clinical=is_mane_plus_clinical,
            mane_status=mane_status,
            window_name=window_name,
            window_start_offset=window_start_offset,
            window_end_offset=window_end_offset,
            promoter_upstream=abs(window_start_offset),
            promoter_downstream=abs(window_end_offset),
            promoter_start=promoter_start,
            promoter_end=promoter_end,
            tf_target=tf_target.upper(),
            biosample=biosample,
            encode_experiment_accession=encode_experiment_accession,
            encode_file_accession=encode_file_accession,
            treatment_context=treatment_context,
            genetic_modification_context=genetic_modification_context,
            processing_tier=processing_tier,
            is_authenticated=False,
            auth_status="SOURCE_UNAVAILABLE",
            observation_status="SOURCE_UNAVAILABLE",
            evidence_status=EvidenceInterpretation.SOURCE_UNAVAILABLE,
            exclusion_reason=f"Primary BED file {encode_file_accession} not found on local disk"
        )

    # 3. Infer or verify chromosome
    chrom_to_search = chromosome
    if not chrom_to_search:
        gene_chrom_map = {
            "CDK4": "chr12", "TFF1": "chr21", "TP53": "chr17",
            "CDKN1A": "chr6", "CTCF": "chr16", "MYC": "chr8", "CCND1": "chr11"
        }
        chrom_to_search = gene_chrom_map.get(gene_symbol.upper(), "")

    # 3b. Resolve transcript version against authoritative pinned references
    target_transcript_id = transcript_id
    resolution_note = None
    if transcript_id:
        resolved_tx, res_note = resolve_transcript_version(gene_symbol, transcript_id)
        if resolved_tx:
            target_transcript_id = resolved_tx
            resolution_note = res_note

    # 4. Strict 12-point Source-Identity and Data Provenance Contract
    contract_outcome = verify_source_identity_contract(
        file_accession=encode_file_accession,
        bed_path=bed_file,
        expected_experiment=encode_experiment_accession,
        expected_tf=tf_target,
        expected_biosample=biosample,
        expected_assembly=assembly,
        gene_symbol=gene_symbol,
        transcript_id=target_transcript_id,
        chromosome=chrom_to_search,
        strand=strand,
        canonical_tss=canonical_tss,
        expected_treatment=treatment_context if explicit_treatment else None,
        expected_modification=genetic_modification_context if explicit_modification else None
    )

    if not contract_outcome.is_fully_authenticated:
        if contract_outcome.local_parsing in [LocalParsingStatus.CORRUPTED_OR_TRUNCATED, LocalParsingStatus.UNREADABLE]:
            auth_st = "FAIL_CLOSED_CORRUPTED_CACHE"
        elif contract_outcome.local_parsing in [LocalParsingStatus.FILE_NOT_FOUND]:
            auth_st = "SOURCE_UNAVAILABLE"
        elif contract_outcome.authoritative_identity == AuthoritativeIdentityStatus.CHECKSUM_MISMATCH:
            auth_st = "CHECKSUM_MISMATCH"
        elif contract_outcome.authoritative_identity == AuthoritativeIdentityStatus.METADATA_UNAVAILABLE:
            auth_st = "METADATA_UNAVAILABLE"
        elif contract_outcome.authoritative_identity == AuthoritativeIdentityStatus.STATUS_REVOKED_OR_ARCHIVED:
            auth_st = "STATUS_REVOKED_OR_ARCHIVED"
        elif contract_outcome.biological_relationship != BiologicalRelationshipStatus.RELATIONSHIP_VALID:
            auth_st = contract_outcome.biological_relationship.value
        else:
            auth_st = "SOURCE_NOT_VERIFIED"

        return EvidenceConfiguration(
            config_id=config_id,
            gene_symbol=gene_symbol.upper(),
            transcript_id=target_transcript_id,
            canonical_tss=canonical_tss,
            strand=strand,
            assembly=assembly,
            chromosome=chrom_to_search,
            is_mane_select=is_mane_select,
            is_mane_plus_clinical=is_mane_plus_clinical,
            mane_status=mane_status,
            window_name=window_name,
            window_start_offset=window_start_offset,
            window_end_offset=window_end_offset,
            promoter_upstream=abs(window_start_offset),
            promoter_downstream=abs(window_end_offset),
            promoter_start=promoter_start,
            promoter_end=promoter_end,
            tf_target=tf_target.upper(),
            biosample=biosample,
            encode_experiment_accession=encode_experiment_accession,
            encode_file_accession=encode_file_accession,
            treatment_context=treatment_context,
            genetic_modification_context=genetic_modification_context,
            processing_tier=processing_tier,
            is_authenticated=False,
            auth_status=auth_st,
            observation_status="SOURCE_NOT_VERIFIED",
            evidence_status=EvidenceInterpretation.SOURCE_NOT_VERIFIED,
            exclusion_reason=contract_outcome.exclusion_reason
        )

    file_sha256 = contract_outcome.actual_sha256
    auth_treat = contract_outcome.token.treatment_context if contract_outcome.token else treatment_context
    auth_mod = contract_outcome.token.genetic_modification_context if contract_outcome.token else genetic_modification_context

    if not chrom_to_search and contract_outcome.token and contract_outcome.token.chromosome:
        chrom_to_search = contract_outcome.token.chromosome

    # 5. Scan the authentic primary BED file across the promoter interval
    overlapping_peaks = scan_bed_for_promoter_overlap(
        bed_file,
        chrom_to_search,
        promoter_start,
        promoter_end
    )

    if overlapping_peaks:
        top_peak = overlapping_peaks[0]
        return EvidenceConfiguration(
            config_id=config_id,
            gene_symbol=gene_symbol.upper(),
            transcript_id=target_transcript_id,
            canonical_tss=canonical_tss,
            strand=strand,
            assembly=assembly,
            chromosome=chrom_to_search,
            is_mane_select=is_mane_select,
            is_mane_plus_clinical=is_mane_plus_clinical,
            mane_status=mane_status,
            window_name=window_name,
            window_start_offset=window_start_offset,
            window_end_offset=window_end_offset,
            promoter_upstream=abs(window_start_offset),
            promoter_downstream=abs(window_end_offset),
            promoter_start=promoter_start,
            promoter_end=promoter_end,
            tf_target=tf_target.upper(),
            biosample=biosample,
            encode_experiment_accession=encode_experiment_accession,
            encode_file_accession=encode_file_accession,
            treatment_context=auth_treat,
            genetic_modification_context=auth_mod,
            processing_tier=processing_tier,
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            source_file_sha256=file_sha256,
            provenance_token=contract_outcome.token,
            is_synthetic=False,
            peak_chromosome=top_peak["chromosome"],
            source_peak_start=top_peak["start"],
            source_peak_end=top_peak["end"],
            source_peak_signal=top_peak["signal"],
            source_peak_summit=top_peak["summit"],
            all_overlapping_peaks=overlapping_peaks,
            has_overlap=True,
            overlap_bp=top_peak["overlap_bp"],
            overlap_percentage=round(top_peak["overlap_percentage"], 2),
            observation_status="VERIFIED_POSITIVE_OBSERVATION",
            evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
            exclusion_reason=None
        )
    else:
        # File scanned completely; no peak detected in promoter interval
        return EvidenceConfiguration(
            config_id=config_id,
            gene_symbol=gene_symbol.upper(),
            transcript_id=target_transcript_id,
            canonical_tss=canonical_tss,
            strand=strand,
            assembly=assembly,
            chromosome=chrom_to_search,
            is_mane_select=is_mane_select,
            is_mane_plus_clinical=is_mane_plus_clinical,
            mane_status=mane_status,
            window_name=window_name,
            window_start_offset=window_start_offset,
            window_end_offset=window_end_offset,
            promoter_upstream=abs(window_start_offset),
            promoter_downstream=abs(window_end_offset),
            promoter_start=promoter_start,
            promoter_end=promoter_end,
            tf_target=tf_target.upper(),
            biosample=biosample,
            encode_experiment_accession=encode_experiment_accession,
            encode_file_accession=encode_file_accession,
            treatment_context=auth_treat,
            genetic_modification_context=auth_mod,
            processing_tier=processing_tier,
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            source_file_sha256=file_sha256,
            provenance_token=contract_outcome.token,
            is_synthetic=False,
            has_overlap=False,
            overlap_bp=0,
            overlap_percentage=0.0,
            observation_status="VERIFIED_NO_PEAK_DETECTED",
            evidence_status=EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED,
            exclusion_reason=None
        )


def create_synthetic_configuration_for_testing(
    config_id: str,
    gene_symbol: str,
    transcript_id: str,
    canonical_tss: int,
    strand: str,
    window_name: str = "Standard",
    window_start_offset: int = -1000,
    window_end_offset: int = 200,
    tf_target: str = "E2F1",
    biosample: str = "K562",
    encode_experiment_accession: str = "ENCSR000EXP",
    encode_file_accession: str = "ENCFF000FIL",
    treatment_context: str = "Untreated",
    genetic_modification_context: str = "None",
    processing_tier: str = "conservative_idr_peaks",
    is_authenticated: bool = True,
    auth_status: str = "AUTHENTICATED_PRIMARY_BED",
    has_overlap: bool = True,
    overlap_bp: int = 100,
    overlap_percentage: float = 50.0,
    evidence_status: EvidenceInterpretation = EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
    exclusion_reason: Optional[str] = None
) -> EvidenceConfiguration:
    """
    Explicitly labeled constructor for synthetic unit testing only.
    MUST NOT be used in researcher-facing production workflows.
    """
    promoter_start, promoter_end = calculate_signed_promoter_coordinates(
        canonical_tss,
        strand,
        window_start_offset,
        window_end_offset,
        is_bed0=True
    )
    if has_overlap:
        source_peak_start = promoter_start + 10
        source_peak_end = promoter_start + 10 + max(1, overlap_bp)
    else:
        source_peak_start = None
        source_peak_end = None

    syn_token = create_synthetic_provenance_token_for_testing(
        file_accession=encode_file_accession,
        experiment_accession=encode_experiment_accession,
        target_tf=tf_target,
        biosample=biosample,
        synthetic_note="Synthetic unit test harness"
    )

    return EvidenceConfiguration(
        config_id=config_id,
        gene_symbol=gene_symbol.upper(),
        transcript_id=transcript_id,
        canonical_tss=canonical_tss,
        strand=strand,
        window_name=window_name,
        window_start_offset=window_start_offset,
        window_end_offset=window_end_offset,
        promoter_upstream=abs(window_start_offset),
        promoter_downstream=abs(window_end_offset),
        promoter_start=promoter_start,
        promoter_end=promoter_end,
        tf_target=tf_target.upper(),
        biosample=biosample,
        encode_experiment_accession=encode_experiment_accession,
        encode_file_accession=encode_file_accession,
        treatment_context=treatment_context,
        genetic_modification_context=genetic_modification_context,
        processing_tier=processing_tier,
        is_authenticated=is_authenticated,
        auth_status=auth_status,
        provenance_token=syn_token,
        is_synthetic=True,
        synthetic_note="Synthetic unit test harness",
        source_peak_start=source_peak_start,
        source_peak_end=source_peak_end,
        has_overlap=has_overlap,
        overlap_bp=overlap_bp,
        overlap_percentage=overlap_percentage,
        observation_status=evidence_status.value,
        evidence_status=evidence_status,
        exclusion_reason=exclusion_reason,
        notes=["SYNTHETIC_UNIT_TEST_RECORD"]
    )
