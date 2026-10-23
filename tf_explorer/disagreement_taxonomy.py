"""
Disagreement Taxonomy and Core Data Models for TF-Explorer Phase 4B-R.

Defines standardized, reproducible classifications explaining disagreement
in promoter-associated transcription-factor evidence across transcripts,
experiments, biological contexts, and analysis settings.
"""

from enum import Enum
from pathlib import Path
from typing import Dict, Any, List, Optional, Set, Tuple
from dataclasses import dataclass, field, asdict

from tf_explorer.disagreement_provenance import (
    scan_bed_for_promoter_overlap,
    normalize_tf_name,
    normalize_biosample_name,
    compare_treatment_contexts,
    compare_genetic_modifications
)


class DisagreementClassification(str, Enum):
    """
    Standardized, reproducible cause-of-disagreement classifications.
    Non-mutually exclusive: multiple interacting causes may be assigned.
    """
    TRANSCRIPT_GEOMETRY_DEPENDENT = "TRANSCRIPT_GEOMETRY_DEPENDENT"
    PROMOTER_WINDOW_DEPENDENT = "PROMOTER_WINDOW_DEPENDENT"
    EXPERIMENT_DISCORDANT = "EXPERIMENT_DISCORDANT"
    BIOLOGICAL_CONTEXT_DIFFERENT = "BIOLOGICAL_CONTEXT_DIFFERENT"
    PEAK_PROCESSING_DEPENDENT = "PEAK_PROCESSING_DEPENDENT"
    INSUFFICIENT_COMPARABLE_DATA = "INSUFFICIENT_COMPARABLE_DATA"
    STABLE_WITHIN_EVALUATED_CHOICES = "STABLE_WITHIN_EVALUATED_CHOICES"


class EvidenceInterpretation(str, Enum):
    """
    Evidence observation and interpretation status for an individual configuration.
    Reflects verifiable physical peak detection or verified absence in primary data.
    """
    VERIFIED_POSITIVE_OBSERVATION = "VERIFIED_POSITIVE_OBSERVATION"
    VERIFIED_NO_PEAK_DETECTED = "VERIFIED_NO_PEAK_DETECTED"
    SOURCE_NOT_VERIFIED = "SOURCE_NOT_VERIFIED"
    SOURCE_UNAVAILABLE = "SOURCE_UNAVAILABLE"
    CONTEXT_INCOMPATIBLE = "CONTEXT_INCOMPATIBLE"
    INSUFFICIENT_COMPARABLE_DATA = "INSUFFICIENT_COMPARABLE_DATA"

    # Backward compatibility aliases
    POSITIVE_EVIDENCE = "VERIFIED_POSITIVE_OBSERVATION"
    NEGATIVE_EVIDENCE = "VERIFIED_NO_PEAK_DETECTED"
    EXCLUDED_UNAUTHENTICATED = "SOURCE_NOT_VERIFIED"
    EXCLUDED_INSUFFICIENT_DATA = "INSUFFICIENT_COMPARABLE_DATA"


class DimensionType(str, Enum):
    """Dimensions of variation across regulatory evidence configurations."""
    TRANSCRIPT_TSS = "TRANSCRIPT_TSS"
    PROMOTER_WINDOW = "PROMOTER_WINDOW"
    EXPERIMENT_ASSAY = "EXPERIMENT_ASSAY"
    BIOLOGICAL_CONTEXT = "BIOLOGICAL_CONTEXT"
    TREATMENT_PERTURBATION = "TREATMENT_PERTURBATION"
    GENETIC_MODIFICATION = "GENETIC_MODIFICATION"
    PEAK_PROCESSING = "PEAK_PROCESSING"
    DATA_AVAILABILITY = "DATA_AVAILABILITY"


def calculate_signed_promoter_coordinates(
    tss: int,
    strand: str,
    window_start_offset: int = -1000,
    window_end_offset: int = 200,
    is_bed0: bool = True
) -> Tuple[int, int]:
    """
    Calculates genomic promoter window [promoter_start, promoter_end) in BED0 coordinates
    for arbitrary signed offsets relative to TSS.
    
    Offsets are in 5' -> 3' transcription coordinates:
      offset < 0 is upstream (5') of TSS
      offset = 0 is at TSS
      offset > 0 is downstream (3') of TSS into the gene body
      
    window_start_offset must be strictly less than window_end_offset.
    
    Examples:
      Standard core: [-500, +100]
      Downstream only: [0, +500]
      Downstream shifted: [+50, +500]
      Upstream only: [-1000, -200]
    """
    if window_start_offset >= window_end_offset:
        raise ValueError(
            f"window_start_offset ({window_start_offset}) must be strictly less than "
            f"window_end_offset ({window_end_offset})"
        )

    strand = strand.strip()
    if strand == "+":
        # Positive strand: transcription proceeds in increasing genomic coordinates.
        # 1-based TSS T corresponds to 0-based boundary T - 1 if is_bed0.
        t_base = tss - 1 if is_bed0 else tss
        p_start = t_base + window_start_offset
        p_end = t_base + window_end_offset
    elif strand == "-":
        # Negative strand: transcription proceeds in decreasing genomic coordinates.
        # 1-based TSS T corresponds to 0-based right-open boundary T.
        # Upstream (window_start_offset < 0) corresponds to higher genomic coordinates: T - window_start_offset
        # Downstream (window_end_offset > 0) corresponds to lower genomic coordinates: T - window_end_offset
        p_start = tss - window_end_offset
        p_end = tss - window_start_offset
    else:
        raise ValueError(f"Invalid strand: '{strand}'. Must be '+' or '-'.")

    return p_start, p_end


@dataclass
class EvidenceConfiguration:
    """
    Evaluated configuration in the Multidimensional Evidence Matrix.
    Represents one coordinate: Transcript x Promoter Window x Experiment x Biological Context.
    Derived deterministically from source genomic and primary BED evidence.
    """
    config_id: str
    gene_symbol: str
    ncbi_gene_id: Optional[str] = None
    ensembl_gene_id: Optional[str] = None
    
    # Transcript dimensions
    transcript_id: str = ""
    transcript_version: Optional[str] = None
    canonical_tss: int = 0
    strand: str = "+"
    assembly: str = "GRCh38"
    chromosome: str = ""
    is_mane_select: bool = False
    is_mane_plus_clinical: bool = False
    mane_status: str = "MANE_UNVERIFIED"
    
    # Promoter window dimensions
    window_name: str = "Standard"
    window_start_offset: int = -1000
    window_end_offset: int = 200
    promoter_upstream: int = 1000
    promoter_downstream: int = 200
    promoter_start: int = 0
    promoter_end: int = 0
    
    # Experimental & Assay dimensions
    tf_target: str = ""
    biosample: str = ""
    encode_experiment_accession: str = ""
    encode_file_accession: str = ""
    treatment_context: str = "Untreated"
    genetic_modification_context: str = "None"
    processing_tier: str = "conservative_idr_peaks"
    
    # Source Peak & Authentication
    is_authenticated: bool = False
    auth_status: str = "SOURCE_NOT_VERIFIED"
    source_file_sha256: Optional[str] = None
    peak_chromosome: Optional[str] = None
    source_peak_start: Optional[int] = None
    source_peak_end: Optional[int] = None
    source_peak_signal: Optional[float] = None
    source_peak_summit: Optional[int] = None
    all_overlapping_peaks: List[Dict[str, Any]] = field(default_factory=list)
    provenance_token: Optional[Any] = None
    is_synthetic: bool = False
    synthetic_note: Optional[str] = None
    
    # Overlap arithmetic (BED0 half-open)
    has_overlap: bool = False
    overlap_bp: int = 0
    overlap_percentage: float = 0.0
    
    # Derived Observation Status
    observation_status: str = "SOURCE_NOT_VERIFIED"
    evidence_status: EvidenceInterpretation = EvidenceInterpretation.SOURCE_NOT_VERIFIED
    exclusion_reason: Optional[str] = None
    notes: List[str] = field(default_factory=list)

    def verify_token_binding(self) -> Tuple[bool, str]:
        """
        Enforces strict cryptographic and biological binding between this EvidenceConfiguration
        and its attached ProvenanceToken.
        Prevents token replay attacks where a valid token from one gene/TF/experiment/file
        is attached to an unrelated configuration.
        """
        if self.provenance_token is None:
            return False, "No provenance token attached to configuration"

        token = self.provenance_token

        # 1. Handle synthetic tokens
        tok_synthetic = bool(getattr(token, "is_synthetic", False))
        if tok_synthetic or self.is_synthetic:
            if not self.is_synthetic:
                return False, "Synthetic provenance token cannot be attached to authentic configuration"
            return True, "Synthetic test configuration bound to synthetic token"

        # 2. Check token internal validity and disk checksums
        valid_tok, tok_reason = token.is_valid_and_unmutated()
        if not valid_tok:
            return False, f"Attached provenance token failed validation: {tok_reason}"

        # 3. Gene identity binding
        if not self.gene_symbol or str(self.gene_symbol).strip() == "":
            return False, "Mandatory biological identity field 'gene_symbol' is missing or empty on configuration"
        cfg_gene = str(self.gene_symbol).strip().upper()
        tok_gene = str(token.gene_symbol).strip().upper()
        if cfg_gene != tok_gene:
            return False, f"Gene symbol mismatch: configuration claims '{cfg_gene}', token was issued for '{tok_gene}'"

        # 4. Exact Versioned Transcript binding
        if not self.transcript_id or str(self.transcript_id).strip() == "":
            return False, "Mandatory biological identity field 'transcript_id' is missing or empty on configuration"
        cfg_tx = str(self.transcript_id).strip()
        tok_tx = str(token.transcript_id).strip()
        if "." not in cfg_tx:
            return False, f"Configuration transcript ID '{cfg_tx}' is unversioned. Verified configurations require exact versioned transcript IDs."
        if cfg_tx != tok_tx:
            return False, f"Transcript ID mismatch: configuration claims '{cfg_tx}', token was issued for '{tok_tx}'"

        # 5. Chromosome binding
        if not self.chromosome or str(self.chromosome).strip() == "":
            return False, "Mandatory biological identity field 'chromosome' is missing or empty on configuration"
        cfg_chrom = str(self.chromosome).strip().lower()
        if not cfg_chrom.startswith("chr"):
            cfg_chrom = f"chr{cfg_chrom}"
        tok_chrom = str(token.chromosome).strip().lower()
        if not tok_chrom.startswith("chr"):
            tok_chrom = f"chr{tok_chrom}"
        if cfg_chrom != tok_chrom:
            return False, f"Chromosome mismatch: configuration claims '{cfg_chrom}', token was issued for '{tok_chrom}'"

        # 6. Strand binding
        if not self.strand or str(self.strand).strip() == "":
            return False, "Mandatory biological identity field 'strand' is missing or empty on configuration"
        if str(self.strand).strip() != str(token.strand).strip():
            return False, f"Strand mismatch: configuration claims '{self.strand}', token was issued for '{token.strand}'"

        # 7. Canonical TSS exact integer binding
        if not self.canonical_tss:
            return False, "Mandatory biological identity field 'canonical_tss' is missing or empty on configuration"
        if self.canonical_tss != token.canonical_tss:
            return False, f"TSS mismatch: configuration claims {self.canonical_tss}, token was issued for {token.canonical_tss}"

        # 8. Genome assembly binding
        if not self.assembly or str(self.assembly).strip() == "":
            return False, "Mandatory biological identity field 'assembly' is missing or empty on configuration"
        if str(self.assembly).strip().upper() != str(token.assembly).strip().upper() or str(self.assembly).strip().upper() != "GRCH38":
            return False, f"Genome assembly mismatch: configuration claims '{self.assembly}', token was issued for '{token.assembly}'"

        # 9. Target TF binding
        if not self.tf_target or str(self.tf_target).strip() == "":
            return False, "Mandatory biological identity field 'tf_target' is missing or empty on configuration"
        cfg_tf = normalize_tf_name(str(self.tf_target))
        tok_tf = normalize_tf_name(str(token.target_tf))
        if cfg_tf != tok_tf:
            return False, f"Target TF mismatch: configuration claims '{self.tf_target}', token was issued for '{token.target_tf}'"

        # 10. Biosample binding
        if not self.biosample or str(self.biosample).strip() == "":
            return False, "Mandatory biological identity field 'biosample' is missing or empty on configuration"
        cfg_bio = normalize_biosample_name(str(self.biosample))
        tok_bio = normalize_biosample_name(str(token.biosample))
        if cfg_bio != tok_bio:
            return False, f"Biosample mismatch: configuration claims '{self.biosample}', token was issued for '{token.biosample}'"

        # 11. ENCODE Experiment accession binding
        if not self.encode_experiment_accession or str(self.encode_experiment_accession).strip() == "":
            return False, "Mandatory biological identity field 'encode_experiment_accession' is missing or empty on configuration"
        cfg_exp = str(self.encode_experiment_accession).strip().upper()
        tok_exp = str(token.experiment_accession).strip().upper()
        if cfg_exp != tok_exp:
            return False, f"Experiment accession mismatch: configuration claims '{cfg_exp}', token was issued for '{tok_exp}'"

        # 12. ENCODE File accession binding
        if not self.encode_file_accession or str(self.encode_file_accession).strip() == "":
            return False, "Mandatory biological identity field 'encode_file_accession' is missing or empty on configuration"
        cfg_file = str(self.encode_file_accession).strip().upper()
        tok_file = str(token.file_accession).strip().upper()
        if cfg_file != tok_file:
            return False, f"File accession mismatch: configuration claims '{cfg_file}', token was issued for '{tok_file}'"

        # 13. Source file SHA-256 binding
        if self.source_file_sha256 is None or str(self.source_file_sha256).strip() == "":
            return False, "Mandatory biological identity field 'source_file_sha256' is missing or empty on configuration"
        cfg_sha = str(self.source_file_sha256).strip().lower()
        tok_sha = str(token.source_file_sha256).strip().lower()
        if cfg_sha != tok_sha:
            return False, f"Source file SHA-256 mismatch: configuration claims '{cfg_sha}', token was issued for '{tok_sha}'"

        # 14. Treatment context binding
        if not self.treatment_context or str(self.treatment_context).strip() == "":
            return False, "Mandatory biological identity field 'treatment_context' is missing or empty on configuration"
        valid_t, msg_t = compare_treatment_contexts(self.treatment_context, token.treatment_context)
        if not valid_t:
            return False, f"Treatment context mismatch: {msg_t}"

        # 15. Genetic modification binding
        if not self.genetic_modification_context or str(self.genetic_modification_context).strip() == "":
            return False, "Mandatory biological identity field 'genetic_modification_context' is missing or empty on configuration"
        valid_m, msg_m = compare_genetic_modifications(self.genetic_modification_context, token.genetic_modification_context)
        if not valid_m:
            return False, f"Genetic modification mismatch: {msg_m}"

        return True, "Configuration strictly bound to provenance token"

    def authenticate_observation(self) -> Tuple[bool, str]:
        """
        Independently authenticates the claimed observation state against the primary BED file bytes.
        Does NOT trust caller-provided has_overlap, overlap_bp, or evidence_status.
        Recomputes signed promoter bounds and scans the authentic BED file.
        """
        if self.is_synthetic:
            return True, "Synthetic test observation"

        if self.provenance_token is None:
            return False, "No provenance token present"

        # Must have valid token binding first
        bound, bind_reason = self.verify_token_binding()
        if not bound:
            return False, f"Cannot authenticate observation: {bind_reason}"

        # 1. Verify promoter coordinates match mathematical signed definition
        try:
            exp_p_start, exp_p_end = calculate_signed_promoter_coordinates(
                self.canonical_tss,
                self.strand,
                self.window_start_offset,
                self.window_end_offset,
                is_bed0=True
            )
        except Exception as e:
            return False, f"Failed to compute promoter coordinates: {e}"

        if self.promoter_start != exp_p_start or self.promoter_end != exp_p_end:
            return False, (
                f"Promoter coordinates mismatch: claimed [{self.promoter_start}, {self.promoter_end}), "
                f"mathematically expected [{exp_p_start}, {exp_p_end}) for TSS={self.canonical_tss}, "
                f"strand={self.strand}, offsets=[{self.window_start_offset}, {self.window_end_offset}]"
            )

        # 2. Check source BED file on disk
        source_path = Path(self.provenance_token.source_file_path)
        if not source_path.exists():
            return False, f"Primary BED file not found on disk at: {source_path}"

        chrom_to_search = self.chromosome or self.provenance_token.chromosome
        if not chrom_to_search:
            return False, "Chromosome not specified for promoter overlap scan"

        # 3. Scan primary BED records directly
        try:
            actual_peaks = scan_bed_for_promoter_overlap(
                source_path,
                chrom_to_search,
                exp_p_start,
                exp_p_end
            )
        except Exception as e:
            return False, f"Failed to scan primary BED records for overlap: {e}"

        expected_has_overlap = bool(actual_peaks)
        if expected_has_overlap:
            expected_status = EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION
            top_peak = actual_peaks[0]
            expected_overlap_bp = top_peak["overlap_bp"]
        else:
            expected_status = EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED
            expected_overlap_bp = 0

        # 4. Check caller-asserted observation fields
        if self.has_overlap != expected_has_overlap:
            return False, (
                f"Observation contradiction: configuration claims has_overlap={self.has_overlap}, "
                f"but primary BED scan derived has_overlap={expected_has_overlap} "
                f"({len(actual_peaks)} overlapping peak(s) found in [{exp_p_start}, {exp_p_end}))"
            )

        if self.evidence_status != expected_status:
            return False, (
                f"Observation contradiction: configuration claims evidence_status='{self.evidence_status.value}', "
                f"but primary BED scan derived evidence_status='{expected_status.value}'"
            )

        if expected_has_overlap:
            if self.overlap_bp != expected_overlap_bp:
                return False, (
                    f"Observation overlap_bp mismatch: configuration claims {self.overlap_bp} bp, "
                    f"primary BED scan derived {expected_overlap_bp} bp"
                )
        else:
            if self.overlap_bp != 0:
                return False, (
                    f"Observation overlap_bp mismatch: configuration claims {self.overlap_bp} bp "
                    f"for a no-peak observation (expected 0 bp)"
                )

        return True, "Observation authenticated against primary BED records"

    def __post_init__(self):
        # Handle Provenance Token strictly fail-closed
        if self.provenance_token is not None:
            valid, reason = self.provenance_token.is_valid_and_unmutated()
            if not valid:
                self.is_authenticated = False
                self.auth_status = "SOURCE_NOT_VERIFIED"
                self.observation_status = "SOURCE_NOT_VERIFIED"
                self.evidence_status = EvidenceInterpretation.SOURCE_NOT_VERIFIED
                self.exclusion_reason = reason
                return

            self.is_synthetic = bool(getattr(self.provenance_token, "is_synthetic", False))
            self.synthetic_note = getattr(self.provenance_token, "synthetic_note", None)

            if self.is_synthetic:
                self.is_authenticated = False
                self.auth_status = "SYNTHETIC_DEMONSTRATION"
                if not self.exclusion_reason:
                    self.exclusion_reason = f"Synthetic/manual configuration excluded from biological consensus ({self.synthetic_note or 'Synthetic unit test harness'})"
                return

            # Strict Token-Configuration Binding
            bound, bind_reason = self.verify_token_binding()
            if not bound:
                self.is_authenticated = False
                self.auth_status = "TOKEN_CONFIGURATION_MISMATCH"
                self.observation_status = "SOURCE_NOT_VERIFIED"
                self.evidence_status = EvidenceInterpretation.SOURCE_NOT_VERIFIED
                self.exclusion_reason = f"Token-configuration binding failed: {bind_reason}"
                return

            # Strict Observation Authentication against primary BED records
            obs_valid, obs_reason = self.authenticate_observation()
            if not obs_valid:
                self.is_authenticated = False
                self.auth_status = "OBSERVATION_CONTRADICTS_PRIMARY_DATA"
                self.observation_status = "SOURCE_NOT_VERIFIED"
                self.evidence_status = EvidenceInterpretation.SOURCE_NOT_VERIFIED
                self.exclusion_reason = f"Observation authentication failed: {obs_reason}"
                return

            self.is_authenticated = True
            self.auth_status = "AUTHENTICATED_PRIMARY_BED"
        else:
            # If no provenance token was provided:
            unverified_statuses = {
                EvidenceInterpretation.SOURCE_UNAVAILABLE,
                EvidenceInterpretation.SOURCE_NOT_VERIFIED,
                EvidenceInterpretation.CONTEXT_INCOMPATIBLE,
                EvidenceInterpretation.INSUFFICIENT_COMPARABLE_DATA
            }
            if self.evidence_status in unverified_statuses:
                self.is_synthetic = False
                self.is_authenticated = False
                if not self.exclusion_reason:
                    self.exclusion_reason = f"Unauthenticated source: {self.auth_status}"
            elif not self.is_authenticated and not self.is_synthetic:
                # Contradictory unauthenticated state: is_authenticated=False and not is_synthetic but claims positive/negative evidence
                self.is_synthetic = False
                self.exclusion_reason = f"Contradictory unauthenticated state (claimed {self.evidence_status.value})"
            else:
                # Manually asserted configuration without sealed provenance token:
                # Strictly isolated as synthetic test mock
                self.is_synthetic = True
                self.is_authenticated = False
                self.auth_status = "SYNTHETIC_DEMONSTRATION"
                self.synthetic_note = "Manually asserted configuration without sealed provenance token"
                self.exclusion_reason = "Synthetic/manual configuration excluded from biological consensus (Manually asserted configuration without sealed provenance token)"

    def has_valid_provenance(self) -> bool:
        """Returns True only if supported by an authentic, unmutated ProvenanceToken."""
        if self.provenance_token is None or self.is_synthetic:
            return False
        valid, _ = self.provenance_token.is_valid_and_unmutated()
        return valid

    def is_verified_biological_observation(self) -> bool:
        """
        Returns True strictly for authentic primary biological observations.
        Excludes synthetic test mocks, unauthenticated configs, missing-token,
        token-mismatched, observation-contradictory, and mutated data.
        """
        if not self.has_valid_provenance() or self.is_synthetic or not self.is_authenticated:
            return False
        if self.auth_status != "AUTHENTICATED_PRIMARY_BED":
            return False
        bound, _ = self.verify_token_binding()
        if not bound:
            return False
        obs_valid, _ = self.authenticate_observation()
        if not obs_valid:
            return False
        return self.evidence_status in [
            EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
            EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED
        ]

    def is_verified_observation(self) -> bool:
        """Returns True if observation is valid within its evaluation domain."""
        if self.is_synthetic:
            return self.evidence_status in [
                EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
                EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED,
                EvidenceInterpretation.POSITIVE_EVIDENCE,
                EvidenceInterpretation.NEGATIVE_EVIDENCE
            ]
        if self.provenance_token is not None:
            valid, _ = self.provenance_token.is_valid_and_unmutated()
            if not valid:
                return False
        return self.is_authenticated and self.evidence_status in [
            EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
            EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED,
            EvidenceInterpretation.POSITIVE_EVIDENCE,
            EvidenceInterpretation.NEGATIVE_EVIDENCE
        ]

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["evidence_status"] = self.evidence_status.value
        if self.provenance_token is not None and hasattr(self.provenance_token, "to_dict"):
            d["provenance_token"] = self.provenance_token.to_dict()
        return d


@dataclass
class DisagreementDiagnosis:
    """
    Structured causal explanation of an identified disagreement between configurations.
    """
    primary_causes: List[DisagreementClassification]
    config_a_id: str
    config_b_id: str
    evidence_a: EvidenceInterpretation
    evidence_b: EvidenceInterpretation
    differing_dimensions: List[DimensionType]
    smallest_observed_difference: str
    scientific_explanation: str
    causal_attribution_caveat: Optional[str] = None
    has_multiple_covarying_dimensions: bool = False
    covarying_dimensions: List[str] = field(default_factory=list)
    has_synthetic_configuration: bool = False

    @property
    def is_multivariate_interaction(self) -> bool:
        """Backward-compatibility alias for has_multiple_covarying_dimensions."""
        return self.has_multiple_covarying_dimensions

    @property
    def interacting_factors(self) -> List[str]:
        """Backward-compatibility alias for covarying_dimensions."""
        return self.covarying_dimensions

    def to_dict(self) -> Dict[str, Any]:
        return {
            "primary_causes": [c.value for c in self.primary_causes],
            "config_a_id": self.config_a_id,
            "config_b_id": self.config_b_id,
            "evidence_a": self.evidence_a.value,
            "evidence_b": self.evidence_b.value,
            "differing_dimensions": [d.value for d in self.differing_dimensions],
            "smallest_observed_difference": self.smallest_observed_difference,
            "scientific_explanation": self.scientific_explanation,
            "causal_attribution_caveat": self.causal_attribution_caveat,
            "has_multiple_covarying_dimensions": self.has_multiple_covarying_dimensions,
            "is_multivariate_interaction": self.has_multiple_covarying_dimensions,
            "covarying_dimensions": self.covarying_dimensions,
            "interacting_factors": self.covarying_dimensions,
            "has_synthetic_configuration": self.has_synthetic_configuration,
        }
