"""
Regulatory Evidence Cause-of-Disagreement Engine for TF-Explorer Phase 4B-R.

Analyzes differences in promoter-associated transcription factor evidence across
transcripts, promoter windows, experiments, biological contexts, and processing settings.
Strictly enforces fail-closed evidence integrity: unauthenticated configurations
cannot contribute positive or negative evidence under any circumstance.
"""

from typing import List, Dict, Any, Optional, Set, Tuple
from dataclasses import dataclass, field

from tf_explorer.disagreement_taxonomy import (
    DisagreementClassification,
    EvidenceInterpretation,
    DimensionType,
    EvidenceConfiguration,
    DisagreementDiagnosis
)
from tf_explorer.disagreement_provenance import (
    compare_treatment_contexts,
    compare_genetic_modifications
)



def classify_peak_profile_concordance(
    config_a: EvidenceConfiguration,
    config_b: EvidenceConfiguration,
    is_both_positive: bool,
    is_both_negative: bool
) -> str:
    """
    Evaluates physical peak profile concordance between two evidence configurations.
    Enforces strict comparability rules:
    - Target TF must match.
    - Assemblies must be GRCh38.
    - Sources must be authenticated.
    - Per-peak coordinate details must be present.
    """
    if not (is_both_positive or is_both_negative):
        return "INCOMPARABLE"
    if is_both_negative:
        return "NO_PEAKS"

    # Biological comparability guards
    if getattr(config_a, "tf_target", None) != getattr(config_b, "tf_target", None):
        return "INCOMPARABLE"
    if getattr(config_a, "assembly", "GRCh38") != "GRCh38" or getattr(config_b, "assembly", "GRCh38") != "GRCh38":
        return "INCOMPARABLE"
    if getattr(config_a, "assembly", "GRCh38") != getattr(config_b, "assembly", "GRCh38"):
        return "INCOMPARABLE"
    # Source authenticity guard: must be authenticated biological data or explicitly synthetic test fixture
    is_syn_a = bool(getattr(config_a, "is_synthetic", False))
    is_syn_b = bool(getattr(config_b, "is_synthetic", False))
    if is_syn_a != is_syn_b:
        return "INCOMPARABLE"
    if not is_syn_a and (not getattr(config_a, "is_authenticated", False) or not getattr(config_b, "is_authenticated", False)):
        return "INCOMPARABLE"

    # Per-peak coordinate extraction
    peaks_a = getattr(config_a, "all_overlapping_peaks", []) or []
    peaks_b = getattr(config_b, "all_overlapping_peaks", []) or []

    if not peaks_a and getattr(config_a, "source_peak_start", None) is not None and getattr(config_a, "source_peak_end", None) is not None:
        peaks_a = [{"chrom": getattr(config_a, "chromosome", "chr"), "start": config_a.source_peak_start, "end": config_a.source_peak_end}]
    if not peaks_b and getattr(config_b, "source_peak_start", None) is not None and getattr(config_b, "source_peak_end", None) is not None:
        peaks_b = [{"chrom": getattr(config_b, "chromosome", "chr"), "start": config_b.source_peak_start, "end": config_b.source_peak_end}]

    # If coordinates are missing, conservative return
    if not peaks_a or not peaks_b:
        return "INCOMPARABLE"

    if len(peaks_a) != len(peaks_b):
        return "DIFFERING_PEAK_COUNTS"

    set_a = {(p.get("chrom"), p.get("start"), p.get("end")) for p in peaks_a}
    set_b = {(p.get("chrom"), p.get("start"), p.get("end")) for p in peaks_b}

    if set_a == set_b:
        return "IDENTICAL_PEAKS"

    shared = set_a & set_b
    if shared:
        return "PARTIALLY_SHARED_PEAKS"

    # Check if intervals overlap each other with boundary offsets
    overlap_count = 0
    for pa in peaks_a:
        ca, sa, ea = pa.get("chrom"), pa.get("start"), pa.get("end")
        for pb in peaks_b:
            cb, sb, eb = pb.get("chrom"), pb.get("start"), pb.get("end")
            if ca == cb and max(sa, sb) < min(ea, eb):
                overlap_count += 1
                break
    if overlap_count == len(peaks_a):
        return "OVERLAPPING_PEAKS_OFFSET"
    elif overlap_count > 0:
        return "PARTIALLY_SHARED_PEAKS"
    else:
        return "DISTINCT_PEAKS"


class DisagreementEngine:
    """
    Evaluates configurations and identifies structured, reproducible causes of disagreement.
    Operates strictly fail-closed: unauthenticated sources are excluded from positive/negative consensus.
    """

    def __init__(self):
        pass

    def evaluate_pair(
        self,
        config_a: EvidenceConfiguration,
        config_b: EvidenceConfiguration,
        allow_synthetic_fixtures: bool = False
    ) -> DisagreementDiagnosis:
        """
        Diagnoses agreement or disagreement between two specific evidence configurations.
        Fails closed if either configuration is unauthenticated, corrupted, or contradictory.
        Separates detection-status concordance from configuration concordance and peak-profile concordance.
        """
        causes: List[DisagreementClassification] = []
        differing_dimensions: List[DimensionType] = []
        covarying_dims: List[str] = []
        caveats: List[str] = []

        valid_auth_statuses = {"AUTHENTICATED_PRIMARY_BED"}
        unverified_evidence_statuses = {
            EvidenceInterpretation.SOURCE_NOT_VERIFIED,
            EvidenceInterpretation.SOURCE_UNAVAILABLE,
            EvidenceInterpretation.CONTEXT_INCOMPATIBLE,
            EvidenceInterpretation.INSUFFICIENT_COMPARABLE_DATA
        }

        # 1. Enforce Fail-Closed Source Authentication & Detect Contradictory / Mutated States
        has_provenance_a = config_a.has_valid_provenance() if hasattr(config_a, "has_valid_provenance") else False
        has_provenance_b = config_b.has_valid_provenance() if hasattr(config_b, "has_valid_provenance") else False
        is_synthetic_a = bool(getattr(config_a, "is_synthetic", False))
        is_synthetic_b = bool(getattr(config_b, "is_synthetic", False))
        is_syn_a = is_synthetic_a or not has_provenance_a
        is_syn_b = is_synthetic_b or not has_provenance_b
        is_synthetic_pair = is_syn_a or is_syn_b

        # Detect contradictory unauthenticated states (unauthenticated and not synthetic, but asserting positive/negative evidence)
        is_contra_a = (not config_a.is_authenticated and not is_synthetic_a) and (config_a.evidence_status not in unverified_evidence_statuses)
        is_contra_b = (not config_b.is_authenticated and not is_synthetic_b) and (config_b.evidence_status not in unverified_evidence_statuses)

        if is_contra_a or is_contra_b:
            causes.append(DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA)
            differing_dimensions.append(DimensionType.DATA_AVAILABILITY)
            missing_info = []
            if is_contra_a:
                missing_info.append(f"{config_a.config_id} (Contradictory unauthenticated state (claimed {config_a.evidence_status.value}))")
            if is_contra_b:
                missing_info.append(f"{config_b.config_id} (Contradictory unauthenticated state (claimed {config_b.evidence_status.value}))")
            smallest_diff = f"Evidence unverified, mutated, or missing in: {'; '.join(missing_info)}"
            return DisagreementDiagnosis(
                primary_causes=causes,
                config_a_id=config_a.config_id,
                config_b_id=config_b.config_id,
                evidence_a=config_a.evidence_status,
                evidence_b=config_b.evidence_status,
                differing_dimensions=differing_dimensions,
                smallest_observed_difference=smallest_diff,
                scientific_explanation="Comparison cannot establish biological or geometric concordance because one or both configurations have contradictory or unverified source status.",
                causal_attribution_caveat="Unverified, synthetic, or absent source data must never be interpreted as positive or negative evidence.",
                has_synthetic_configuration=is_synthetic_a or is_synthetic_b,
                detection_status_concordance="UNVERIFIED",
                configuration_concordance=False,
                peak_profile_concordance="INCOMPARABLE",
                interpretive_uncertainty="Unauthenticated configuration asserting unverified observation status."
            )

        # Detect genuine unverified / unavailable source data
        is_unauth_a = ((not config_a.is_authenticated and not is_synthetic_a) or (config_a.evidence_status in unverified_evidence_statuses))
        is_unauth_b = ((not config_b.is_authenticated and not is_synthetic_b) or (config_b.evidence_status in unverified_evidence_statuses))

        if not is_unauth_a and hasattr(config_a, "verify_token_binding") and config_a.provenance_token is not None and not getattr(config_a, "is_synthetic", False):
            bound_a, bound_reason_a = config_a.verify_token_binding()
            if not bound_a:
                is_unauth_a = True
                config_a.exclusion_reason = f"Token-configuration binding failed: {bound_reason_a}"

        if not is_unauth_b and hasattr(config_b, "verify_token_binding") and config_b.provenance_token is not None and not getattr(config_b, "is_synthetic", False):
            bound_b, bound_reason_b = config_b.verify_token_binding()
            if not bound_b:
                is_unauth_b = True
                config_b.exclusion_reason = f"Token-configuration binding failed: {bound_reason_b}"

        if not is_unauth_a and hasattr(config_a, "authenticate_observation") and config_a.provenance_token is not None and not getattr(config_a, "is_synthetic", False):
            obs_a_ok, obs_reason_a = config_a.authenticate_observation()
            if not obs_a_ok:
                is_unauth_a = True
                config_a.exclusion_reason = f"Observation authentication failed: {obs_reason_a}"

        if not is_unauth_b and hasattr(config_b, "authenticate_observation") and config_b.provenance_token is not None and not getattr(config_b, "is_synthetic", False):
            obs_b_ok, obs_reason_b = config_b.authenticate_observation()
            if not obs_b_ok:
                is_unauth_b = True
                config_b.exclusion_reason = f"Observation authentication failed: {obs_reason_b}"

        if is_unauth_a or is_unauth_b:
            causes.append(DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA)
            differing_dimensions.append(DimensionType.DATA_AVAILABILITY)
            missing_info = []
            if is_unauth_a:
                missing_info.append(f"{config_a.config_id} ({config_a.exclusion_reason or config_a.auth_status})")
            if is_unauth_b:
                missing_info.append(f"{config_b.config_id} ({config_b.exclusion_reason or config_b.auth_status})")
            smallest_diff = f"Evidence unverified, mutated, synthetic, or missing in: {'; '.join(missing_info)}"
            return DisagreementDiagnosis(
                primary_causes=causes,
                config_a_id=config_a.config_id,
                config_b_id=config_b.config_id,
                evidence_a=config_a.evidence_status,
                evidence_b=config_b.evidence_status,
                differing_dimensions=differing_dimensions,
                smallest_observed_difference=smallest_diff,
                scientific_explanation=(
                    "Comparison cannot establish biological or geometric concordance because one or both "
                    "configurations lack authenticated, unmutated primary data with valid provenance tokens. "
                    "TF-Explorer fails closed on unverified sources."
                ),
                causal_attribution_caveat="Unverified, synthetic, or absent source data must never be interpreted as positive or negative evidence.",
                has_synthetic_configuration=getattr(config_a, "is_synthetic", False) or getattr(config_b, "is_synthetic", False),
                detection_status_concordance="UNVERIFIED",
                configuration_concordance=False,
                peak_profile_concordance="INCOMPARABLE",
                interpretive_uncertainty="Unverified, missing, or tampered source data."
            )

        # 2. Asymmetric Provenance Check: Authentic Biological Observation vs Synthetic/Unverified
        # An authenticated observation can NEVER be compared against a synthetic mock or unverified configuration.
        if has_provenance_a != has_provenance_b:
            causes.append(DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA)
            differing_dimensions.append(DimensionType.DATA_AVAILABILITY)
            missing_info = []
            if not has_provenance_a:
                missing_info.append(f"{config_a.config_id} (Unverified/synthetic configuration cannot be compared with authentic evidence)")
            if not has_provenance_b:
                missing_info.append(f"{config_b.config_id} (Unverified/synthetic configuration cannot be compared with authentic evidence)")
            smallest_diff = f"Evidence unverified, synthetic, or missing in: {'; '.join(missing_info)}"
            return DisagreementDiagnosis(
                primary_causes=causes,
                config_a_id=config_a.config_id,
                config_b_id=config_b.config_id,
                evidence_a=config_a.evidence_status,
                evidence_b=config_b.evidence_status,
                differing_dimensions=differing_dimensions,
                smallest_observed_difference=smallest_diff,
                scientific_explanation="Comparison cannot establish biological concordance or discordance because authenticated observations cannot be compared against unverified or synthetic test data.",
                causal_attribution_caveat="Authenticated biological evidence cannot be combined with synthetic mock data.",
                has_synthetic_configuration=True,
                detection_status_concordance="UNVERIFIED",
                configuration_concordance=False,
                peak_profile_concordance="INCOMPARABLE",
                interpretive_uncertainty="Asymmetric provenance: comparison of authenticated observation against unverified/synthetic data."
            )

        # 3. Detect Internal Contradictions on Identical Sources/Parameters
        same_exp = bool(
            config_a.encode_experiment_accession
            and config_a.encode_experiment_accession == config_b.encode_experiment_accession
        )
        same_file = (config_a.encode_file_accession == config_b.encode_file_accession)
        same_gene = (config_a.gene_symbol.strip().upper() == config_b.gene_symbol.strip().upper() if config_a.gene_symbol and config_b.gene_symbol else True)
        same_tss = (config_a.canonical_tss == config_b.canonical_tss and config_a.strand == config_b.strand)
        same_win = (
            config_a.window_start_offset == config_b.window_start_offset
            and config_a.window_end_offset == config_b.window_end_offset
            and config_a.promoter_upstream == config_b.promoter_upstream
            and config_a.promoter_downstream == config_b.promoter_downstream
        )
        same_bio = (
            config_a.biosample.strip().lower() == config_b.biosample.strip().lower()
            if config_a.biosample and config_b.biosample else True
        )
        treat_match, _ = compare_treatment_contexts(config_a.treatment_context, config_b.treatment_context)
        mod_match, _ = compare_genetic_modifications(config_a.genetic_modification_context, config_b.genetic_modification_context)
        status_discordant = (config_a.evidence_status != config_b.evidence_status)

        # If identical parameters assert opposite evidence statuses -> Contradictory Data Integrity Inconsistency!
        if same_exp and same_file and same_gene and same_tss and same_win and same_bio and treat_match and mod_match and status_discordant:
            causes.append(DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA)
            differing_dimensions.append(DimensionType.DATA_AVAILABILITY)
            return DisagreementDiagnosis(
                primary_causes=causes,
                config_a_id=config_a.config_id,
                config_b_id=config_b.config_id,
                evidence_a=config_a.evidence_status,
                evidence_b=config_b.evidence_status,
                differing_dimensions=differing_dimensions,
                smallest_observed_difference="Contradictory evidence statuses asserted on identical experimental accession, file, and coordinates.",
                scientific_explanation="Comparison cannot proceed because identical experimental and coordinate parameters assert mutually contradictory observation statuses. This indicates an internal data-integrity inconsistency or corrupted record.",
                causal_attribution_caveat="Primary data cannot support conflicting binding observations for identical inputs.",
                has_synthetic_configuration=is_synthetic_pair,
                detection_status_concordance="DISCORDANT",
                configuration_concordance=True,
                peak_profile_concordance="INCOMPARABLE",
                interpretive_uncertainty="Internal data-integrity inconsistency: conflicting observation statuses for identical source and coordinate inputs."
            )

        # 4. Check for Synthetic / Manual Configurations
        # Synthetic configurations claiming concordant evidence cannot establish biological consensus.
        if is_synthetic_pair and (config_a.evidence_status == config_b.evidence_status):
            causes.append(DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA)
            differing_dimensions.append(DimensionType.DATA_AVAILABILITY)
            missing_info = []
            if is_syn_a:
                missing_info.append(f"{config_a.config_id} (Synthetic/manual configuration excluded from biological consensus ({getattr(config_a, 'synthetic_note', None) or 'synthetic mock'}))")
            if is_syn_b:
                missing_info.append(f"{config_b.config_id} (Synthetic/manual configuration excluded from biological consensus ({getattr(config_b, 'synthetic_note', None) or 'synthetic mock'}))")
            smallest_diff = f"Evidence unverified, synthetic, or missing in: {'; '.join(missing_info)}"
            return DisagreementDiagnosis(
                primary_causes=causes,
                config_a_id=config_a.config_id,
                config_b_id=config_b.config_id,
                evidence_a=config_a.evidence_status,
                evidence_b=config_b.evidence_status,
                differing_dimensions=differing_dimensions,
                smallest_observed_difference=smallest_diff,
                scientific_explanation="Synthetic or manually asserted test configurations cannot establish biological concordance. TF-Explorer fails closed on non-biological evidence.",
                causal_attribution_caveat="Synthetic test configurations must never be interpreted as verified biological evidence.",
                has_synthetic_configuration=True,
                detection_status_concordance="UNVERIFIED",
                configuration_concordance=False,
                peak_profile_concordance="INCOMPARABLE",
                interpretive_uncertainty="Synthetic configurations cannot establish biological consensus."
            )

        # 5. Multi-Attribute Dimensional Dissection
        # Dimension A: Transcript TSS Geometry
        tss_diff = abs(config_a.canonical_tss - config_b.canonical_tss)
        transcript_differs = tss_diff > 0 or (config_a.transcript_id != config_b.transcript_id)
        if transcript_differs:
            differing_dimensions.append(DimensionType.TRANSCRIPT_TSS)
            covarying_dims.append(f"TSS displacement: {tss_diff:,} bp ({config_a.transcript_id} vs {config_b.transcript_id})")

        # Dimension B: Promoter Window Boundaries
        window_differs = not same_win
        if window_differs:
            differing_dimensions.append(DimensionType.PROMOTER_WINDOW)
            win_desc_a = f"[{config_a.window_start_offset}, {config_a.window_end_offset}]"
            win_desc_b = f"[{config_b.window_start_offset}, {config_b.window_end_offset}]"
            covarying_dims.append(f"Window definition: {win_desc_a} vs {win_desc_b}")

        # Dimension C: Peak Processing Tiers (Same experiment, different processing files)
        file_differs = config_a.encode_file_accession != config_b.encode_file_accession
        processing_differs = same_exp and file_differs
        if processing_differs:
            differing_dimensions.append(DimensionType.PEAK_PROCESSING)
            covarying_dims.append(f"Processing tier/file: {config_a.encode_file_accession} vs {config_b.encode_file_accession}")

        # Dimension D: Biological Context (Biosample, Treatment, Genetic Modification)
        bio_differs = not same_bio
        treat_differs = not treat_match
        mod_differs = not mod_match
        context_differs = bio_differs or treat_differs or mod_differs
        if context_differs:
            differing_dimensions.append(DimensionType.BIOLOGICAL_CONTEXT)
            diff_ctx = []
            if bio_differs:
                diff_ctx.append(f"biosample: {config_a.biosample} vs {config_b.biosample}")
            if treat_differs:
                diff_ctx.append(f"treatment: {config_a.treatment_context} vs {config_b.treatment_context}")
            if mod_differs:
                diff_ctx.append(f"genetic_modification: {config_a.genetic_modification_context} vs {config_b.genetic_modification_context}")
            covarying_dims.append("; ".join(diff_ctx))

        # Dimension E: Independent Experiment Assay
        exp_differs = bool(
            config_a.encode_experiment_accession
            and config_b.encode_experiment_accession
            and config_a.encode_experiment_accession != config_b.encode_experiment_accession
        )
        if exp_differs and not context_differs:
            differing_dimensions.append(DimensionType.EXPERIMENT_ASSAY)
            covarying_dims.append(
                f"Independent experiments in identical context: {config_a.encode_experiment_accession} vs {config_b.encode_experiment_accession}"
            )

        # Overall configuration concordance
        configuration_concordance = not (
            transcript_differs or window_differs or processing_differs or context_differs or exp_differs
        )

        # 6. Peak Profile & Detection Status Concordance
        is_both_positive = (
            config_a.evidence_status in (EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION, EvidenceInterpretation.POSITIVE_EVIDENCE)
            and config_b.evidence_status in (EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION, EvidenceInterpretation.POSITIVE_EVIDENCE)
        )
        is_both_negative = (
            config_a.evidence_status in (EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED, EvidenceInterpretation.NEGATIVE_EVIDENCE)
            and config_b.evidence_status in (EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED, EvidenceInterpretation.NEGATIVE_EVIDENCE)
        )
        is_status_concordant = (config_a.evidence_status == config_b.evidence_status) or is_both_positive or is_both_negative
        detection_status_concordance = (
            "BOTH_POSITIVE" if is_both_positive else ("BOTH_NEGATIVE" if is_both_negative else "DISCORDANT")
        )

        # Physical Peak Profile Concordance via rigorous helper
        peak_profile_concordance = classify_peak_profile_concordance(
            config_a, config_b, is_both_positive, is_both_negative
        )

        # 7. Formulate Classification, Explanations, and Caveats
        has_multiple_dims = len(differing_dimensions) > 1
        interpretive_uncertainty: Optional[str] = None

        if is_status_concordant:
            # If all parameters match and peaks match, observation is truly stable across evaluated choices:
            if configuration_concordance and peak_profile_concordance in ("IDENTICAL_PEAKS", "NO_PEAKS", "INCOMPARABLE"):
                causes.append(DisagreementClassification.STABLE_WITHIN_EVALUATED_CHOICES)
                smallest_diff = "None (Concordant evidence status and identical parameters)"
                explanation = f"Both configurations yield concordant {config_a.evidence_status.value}."
                if is_both_negative:
                    caveats.append(
                        "Concordant non-detection indicates no overlapping peak was called in the evaluated window. "
                        "It does not establish biological non-binding across unassayed contexts or transcriptional inactivity."
                    )
            else:
                # Detection status is concordant, BUT parameters or peak profiles differ!
                # We separately record the differing dimensions as the primary causes of difference:
                if transcript_differs:
                    causes.append(DisagreementClassification.TRANSCRIPT_GEOMETRY_DEPENDENT)
                if window_differs:
                    causes.append(DisagreementClassification.PROMOTER_WINDOW_DEPENDENT)
                if processing_differs:
                    causes.append(DisagreementClassification.PEAK_PROCESSING_DEPENDENT)
                if context_differs:
                    causes.append(DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT)
                if exp_differs and not context_differs:
                    causes.append(DisagreementClassification.EXPERIMENT_DISCORDANT)
                
                # If configuration parameters matched, but physical peaks differ:
                if not causes and peak_profile_concordance in ("DISTINCT_PEAKS", "DIFFERING_PEAK_COUNTS", "PARTIALLY_SHARED_PEAKS", "OVERLAPPING_PEAKS_OFFSET"):
                    causes.append(DisagreementClassification.PEAK_PROCESSING_DEPENDENT if not exp_differs else DisagreementClassification.EXPERIMENT_DISCORDANT)

                smallest_diff = covarying_dims[0] if len(covarying_dims) == 1 else (" | ".join(covarying_dims) if covarying_dims else f"Peak profile: {peak_profile_concordance}")
                explanation = (
                    f"Concordant {config_a.evidence_status.value} observed across configurations, "
                    f"but configuration parameters or peak profiles differ in: {smallest_diff}. "
                    f"Peak profile concordance: {peak_profile_concordance}."
                )
                interpretive_uncertainty = (
                    f"Detection status is concordant ({detection_status_concordance}), "
                    f"but analytical parameters or peak profiles vary across evaluated choices."
                )
        else:
            # Discordant evidence statuses (one positive, one negative)
            if transcript_differs:
                causes.append(DisagreementClassification.TRANSCRIPT_GEOMETRY_DEPENDENT)
            if window_differs:
                causes.append(DisagreementClassification.PROMOTER_WINDOW_DEPENDENT)
            if processing_differs:
                causes.append(DisagreementClassification.PEAK_PROCESSING_DEPENDENT)
            if context_differs:
                causes.append(DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT)
            if exp_differs and not context_differs:
                causes.append(DisagreementClassification.EXPERIMENT_DISCORDANT)

            # Fallback logic for discordant statuses without isolated modeled cause:
            if not causes:
                if exp_differs:
                    causes.append(DisagreementClassification.EXPERIMENT_DISCORDANT)
                    smallest_diff = f"Independent experiments in identical context: {config_a.encode_experiment_accession} vs {config_b.encode_experiment_accession}"
                    explanation = "Independent ChIP-seq experiments in the same biological context yield discordant peak detection."
                else:
                    causes.append(DisagreementClassification.UNEXPLAINED_EVIDENCE_DISCORDANCE)
                    smallest_diff = "Unexplained evidence discordance within the same experiment"
                    explanation = "Evidence disagrees across configurations from the same experiment without an isolated coordinate or modeled metadata shift."
                    interpretive_uncertainty = "Discordant observation statuses cannot be explained by standard coordinate or metadata differences."
            else:
                smallest_diff = covarying_dims[0] if len(covarying_dims) == 1 else " | ".join(covarying_dims)
                explanation_parts = []
                for cause in causes:
                    if cause == DisagreementClassification.TRANSCRIPT_GEOMETRY_DEPENDENT:
                        explanation_parts.append(
                            f"Transcript TSS shift of {tss_diff:,} bp alters promoter coordinate boundaries relative to the peak interval."
                        )
                    elif cause == DisagreementClassification.PROMOTER_WINDOW_DEPENDENT:
                        explanation_parts.append(
                            "Promoter window boundary parameter changes exclude or include the peak interval."
                        )
                    elif cause == DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT:
                        explanation_parts.append(
                            "Biological cell line, treatment perturbation, or genetic tag differs between configurations."
                        )
                    elif cause == DisagreementClassification.PEAK_PROCESSING_DEPENDENT:
                        explanation_parts.append(
                            "Peak calling threshold, IDR tier, or replication processing differs within the same experiment."
                        )
                    elif cause == DisagreementClassification.EXPERIMENT_DISCORDANT:
                        explanation_parts.append(
                            "Independent ChIP-seq experiments in the same biological context yield discordant peak detection."
                        )
                explanation = " ".join(explanation_parts)

        # 8. Add Scientific Attribution Caveats
        if is_both_negative:
            caveats.append(
                'Concordant non-detection indicates no overlapping peak was called in the evaluated window. '
                'It does not establish biological non-binding across unassayed contexts or transcriptional inactivity.'
            )
        if DisagreementClassification.PROMOTER_WINDOW_DEPENDENT in causes:
            caveats.append(
                "A change in promoter overlap is an analytical coordinate finding and does not establish differential biological regulation."
            )
        if DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT in causes:
            caveats.append(
                "Differences across cell lines or treatments are observational unless validated with matched isogenic controls."
            )
        if DisagreementClassification.EXPERIMENT_DISCORDANT in causes:
            caveats.append(
                "Disagreement between independent assays in the same context suggests assay sensitivity, antibody variability, or protocol differences."
            )
        if has_multiple_dims:
            caveats.append(
                "MULTIPLE_COVARYING_DIMENSIONS: More than one analytical or biological parameter varies simultaneously. "
                "Causal attribution to a single driver requires holding orthogonal dimensions constant."
            )
        if peak_profile_concordance == "DIFFERING_PEAK_COUNTS":
            caveats.append(
                "PEAK_PROFILE_DIFFERENCE: Overlapping peak counts differ between configurations despite matching detection status."
            )
        elif peak_profile_concordance == "DISTINCT_PEAKS":
            caveats.append(
                "DISTINCT_PEAK_COORDINATES: Overlapping peaks represent distinct physical genomic coordinates."
            )

        if is_synthetic_pair:
            caveats.append(
                "Synthetic computational test fixture: demonstrates parameter sensitivity; not authenticated biological evidence."
            )

        caveat_str = " | ".join(caveats) if caveats else "Evidence reflects evaluated analytical parameter choices."

        return DisagreementDiagnosis(
            primary_causes=causes,
            config_a_id=config_a.config_id,
            config_b_id=config_b.config_id,
            evidence_a=config_a.evidence_status,
            evidence_b=config_b.evidence_status,
            differing_dimensions=differing_dimensions,
            smallest_observed_difference=smallest_diff,
            scientific_explanation=explanation,
            causal_attribution_caveat=caveat_str,
            has_multiple_covarying_dimensions=has_multiple_dims,
            covarying_dimensions=covarying_dims,
            has_synthetic_configuration=is_synthetic_pair,
            detection_status_concordance=detection_status_concordance,
            configuration_concordance=configuration_concordance,
            peak_profile_concordance=peak_profile_concordance,
            interpretive_uncertainty=interpretive_uncertainty
        )

    def evaluate_set(
        self,
        configurations: List[EvidenceConfiguration]
    ) -> Dict[str, Any]:
        """
        Evaluates a set of evidence configurations, performing all pairwise comparisons
        and synthesizing the global consensus and disagreement landscape.
        Distinguishes:
        - Detection-status concordance.
        - Configuration concordance.
        - Peak-profile concordance.
        - Unresolved interpretation.
        - Invalid or unauthenticated evidence.
        """
        if not configurations:
            return {
                "overall_status": "EMPTY_SET",
                "total_configurations_count": 0,
                "positive_configurations_count": 0,
                "negative_configurations_count": 0,
                "excluded_configurations_count": 0,
                "unique_experiments_count": 0,
                "status_concordant": True,
                "detection_status_concordant": True,
                "configuration_concordant": True,
                "peak_profile_concordant": True,
                "has_unresolved_interpretation": False,
                "has_invalid_or_unauthenticated_evidence": False,
                "observed_causes": [],
                "diagnoses": [],
                "scientific_interpretation": "No configurations provided for evaluation."
            }

        valid_configs = [
            c for c in configurations
            if c.is_verified_biological_observation()
        ]
        excluded_configs = [c for c in configurations if c not in valid_configs]

        pos_configs = [
            c for c in valid_configs
            if c.evidence_status == EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION
        ]
        neg_configs = [
            c for c in valid_configs
            if c.evidence_status == EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED
        ]

        unique_experiments = {
            c.encode_experiment_accession for c in valid_configs
            if c.encode_experiment_accession
        }

        # Pairwise evaluations
        diagnoses: List[DisagreementDiagnosis] = []
        observed_causes_set: Set[DisagreementClassification] = set()

        all_pairs_config_concordant = True
        all_pairs_peak_concordant = True
        has_unresolved = False
        has_invalid = len(excluded_configs) > 0

        for i in range(len(configurations)):
            for j in range(i + 1, len(configurations)):
                diag = self.evaluate_pair(configurations[i], configurations[j])

                if DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA in diag.primary_causes:
                    has_invalid = True
                if (DisagreementClassification.UNEXPLAINED_EVIDENCE_DISCORDANCE in diag.primary_causes 
                    or bool(diag.interpretive_uncertainty)):
                    has_unresolved = True

                if not diag.configuration_concordance:
                    all_pairs_config_concordant = False
                if diag.peak_profile_concordance not in ("IDENTICAL_PEAKS", "NO_PEAKS", "INCOMPARABLE"):
                    all_pairs_peak_concordant = False

                # Identify causes other than pure STABLE
                differing_causes = [
                    c for c in diag.primary_causes 
                    if c != DisagreementClassification.STABLE_WITHIN_EVALUATED_CHOICES
                ]

                has_divergence = (
                    len(differing_causes) > 0
                    or not diag.configuration_concordance
                    or diag.peak_profile_concordance in ("DISTINCT_PEAKS", "DIFFERING_PEAK_COUNTS", "PARTIALLY_SHARED_PEAKS", "OVERLAPPING_PEAKS_OFFSET")
                    or diag.detection_status_concordance == "DISCORDANT"
                )

                if has_divergence:
                    diagnoses.append(diag)
                    for c in differing_causes:
                        observed_causes_set.add(c)

        observed_causes = [c.value for c in sorted(observed_causes_set, key=lambda x: x.value)]

        is_detection_concordant = bool(
            (len(pos_configs) > 0 and len(neg_configs) == 0) or
            (len(neg_configs) > 0 and len(pos_configs) == 0)
        )

        if not valid_configs:
            overall_status = "INSUFFICIENT_DATA_FAIL_CLOSED"
            interp = (
                f"All {len(configurations)} evaluated configurations were excluded due to unverified, "
                f"synthetic, or missing primary data. Fail-closed policy prevents assigning positive or negative evidence."
            )
        elif len(pos_configs) > 0 and len(neg_configs) == 0:
            if all_pairs_config_concordant and all_pairs_peak_concordant:
                overall_status = "ROBUST_POSITIVE_EVIDENCE"
                interp = (
                    f"Unanimous positive evidence across {len(pos_configs)} verified configuration(s) "
                    f"spanning {len(unique_experiments)} independent experiment(s) with identical analytical parameters and peak profiles."
                )
            else:
                overall_status = "CONCORDANT_POSITIVE_PARAMETER_SENSITIVE"
                causes_str = ", ".join(observed_causes) if observed_causes else "PEAK_PROFILE_DIFFERENCE"
                interp = (
                    f"Concordant positive evidence observed across {len(pos_configs)} verified configuration(s) "
                    f"spanning {len(unique_experiments)} independent experiment(s), but analytical parameters or peak profiles vary: {causes_str}. "
                    f"Observations must not be assumed to have identical genomic boundaries or peak multiplicity."
                )
        elif len(neg_configs) > 0 and len(pos_configs) == 0:
            if all_pairs_config_concordant:
                overall_status = "ROBUST_NO_PEAK_DETECTED"
                interp = (
                    f"Unanimous absence of overlapping peaks across {len(neg_configs)} verified configuration(s) "
                    f"spanning {len(unique_experiments)} independent experiment(s). "
                    f"Note: Non-detection within the evaluated window does not demonstrate absence of biological regulation or distal enhancer binding."
                )
            else:
                overall_status = "CONCORDANT_NEGATIVE_PARAMETER_SENSITIVE"
                causes_str = ", ".join(observed_causes) if observed_causes else "PARAMETER_DIFFERENCE"
                interp = (
                    f"Concordant non-detection across {len(neg_configs)} verified configuration(s) "
                    f"spanning {len(unique_experiments)} independent experiment(s), evaluated under differing analytical contexts: {causes_str}. "
                    f"Non-detection within evaluated promoter windows does not establish biological inactivity."
                )
        else:
            overall_status = "DISCORDANT_EVIDENCE"
            top_causes = ", ".join(observed_causes) if observed_causes else "MULTIPLE_FACTORS"
            interp = (
                f"Evidence is discordant across evaluated configurations: {len(pos_configs)} verified positive "
                f"vs {len(neg_configs)} verified no peak detected across {len(unique_experiments)} independent experiment(s). "
                f"Observed causes: {top_causes}."
            )

        if excluded_configs:
            interp += f" Note: {len(excluded_configs)} configuration(s) were excluded fail-closed due to unverified or synthetic data."

        return {
            "overall_status": overall_status,
            "total_configurations_count": len(configurations),
            "positive_configurations_count": len(pos_configs),
            "negative_configurations_count": len(neg_configs),
            "excluded_configurations_count": len(excluded_configs),
            "unique_experiments_count": len(unique_experiments),
            "status_concordant": is_detection_concordant,
            "detection_status_concordant": is_detection_concordant,
            "configuration_concordant": bool(all_pairs_config_concordant and len(valid_configs) > 0),
            "peak_profile_concordant": bool(all_pairs_peak_concordant and len(valid_configs) > 0),
            "has_unresolved_interpretation": has_unresolved,
            "has_invalid_or_unauthenticated_evidence": has_invalid,
            "observed_causes": observed_causes,
            "diagnoses": [d.to_dict() for d in diagnoses],
            "scientific_interpretation": interp
        }
