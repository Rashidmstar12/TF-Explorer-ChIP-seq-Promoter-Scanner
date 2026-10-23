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
        config_b: EvidenceConfiguration
    ) -> DisagreementDiagnosis:
        """
        Diagnoses disagreement between two specific evidence configurations.
        Fails closed if either configuration is unauthenticated or has contradictory status.
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
                has_synthetic_configuration=is_synthetic_a or is_synthetic_b
            )

        # Detect genuine unverified / unavailable source data (excluding synthetic test configurations)
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
                has_synthetic_configuration=getattr(config_a, "is_synthetic", False) or getattr(config_b, "is_synthetic", False)
            )

        # Check for synthetic / manual configurations without valid biological provenance
        is_syn_a = getattr(config_a, "is_synthetic", False) or not has_provenance_a
        is_syn_b = getattr(config_b, "is_synthetic", False) or not has_provenance_b

        # Synthetic configurations claiming concordant evidence cannot establish biological consensus
        if (is_syn_a or is_syn_b) and (config_a.evidence_status == config_b.evidence_status):
            causes.append(DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA)
            differing_dimensions.append(DimensionType.DATA_AVAILABILITY)
            missing_info = []
            if is_syn_a:
                missing_info.append(f"{config_a.config_id} (Synthetic/manual configuration excluded from biological consensus ({config_a.synthetic_note or 'synthetic mock'}))")
            if is_syn_b:
                missing_info.append(f"{config_b.config_id} (Synthetic/manual configuration excluded from biological consensus ({config_b.synthetic_note or 'synthetic mock'}))")
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
                has_synthetic_configuration=True
            )

        # 2. Both configurations are authenticated and verified. Check for observation status agreement.
        if config_a.evidence_status == config_b.evidence_status:
            # Full observation concordance (both verified positive or both verified no peak)
            return DisagreementDiagnosis(
                primary_causes=[DisagreementClassification.STABLE_WITHIN_EVALUATED_CHOICES],
                config_a_id=config_a.config_id,
                config_b_id=config_b.config_id,
                evidence_a=config_a.evidence_status,
                evidence_b=config_b.evidence_status,
                differing_dimensions=[],
                smallest_observed_difference="None (Concordant evidence status)",
                scientific_explanation=f"Both configurations yield concordant {config_a.evidence_status.value}.",
                has_synthetic_configuration=is_syn_a or is_syn_b
            )

        # 3. Disagreement exists between authenticated configurations A and B!
        # Dissect each differing dimension:

        # Dimension A: Transcript TSS Geometry
        tss_diff = abs(config_a.canonical_tss - config_b.canonical_tss)
        if tss_diff > 0 or (config_a.transcript_id != config_b.transcript_id):
            differing_dimensions.append(DimensionType.TRANSCRIPT_TSS)
            causes.append(DisagreementClassification.TRANSCRIPT_GEOMETRY_DEPENDENT)
            covarying_dims.append(f"TSS displacement: {tss_diff:,} bp ({config_a.transcript_id} vs {config_b.transcript_id})")

        # Dimension B: Promoter Window Boundaries
        window_differs = (
            config_a.window_start_offset != config_b.window_start_offset
            or config_a.window_end_offset != config_b.window_end_offset
            or config_a.promoter_upstream != config_b.promoter_upstream
            or config_a.promoter_downstream != config_b.promoter_downstream
        )
        if window_differs:
            differing_dimensions.append(DimensionType.PROMOTER_WINDOW)
            causes.append(DisagreementClassification.PROMOTER_WINDOW_DEPENDENT)
            win_desc_a = f"[{config_a.window_start_offset}, {config_a.window_end_offset}]"
            win_desc_b = f"[{config_b.window_start_offset}, {config_b.window_end_offset}]"
            covarying_dims.append(f"Window definition: {win_desc_a} vs {win_desc_b}")

        # Dimension C: Peak Processing Tiers (Same experiment, different processing files)
        same_exp = (
            config_a.encode_experiment_accession
            and config_a.encode_experiment_accession == config_b.encode_experiment_accession
        )
        file_differs = config_a.encode_file_accession != config_b.encode_file_accession

        if same_exp and file_differs:
            differing_dimensions.append(DimensionType.PEAK_PROCESSING)
            causes.append(DisagreementClassification.PEAK_PROCESSING_DEPENDENT)
            covarying_dims.append(f"Processing tier/file: {config_a.encode_file_accession} vs {config_b.encode_file_accession}")

        # Dimension D: Biological Context (Biosample, Treatment, Genetic Modification)
        bio_differs = (
            config_a.biosample.strip().lower() != config_b.biosample.strip().lower()
            if config_a.biosample and config_b.biosample else False
        )
        treat_match, _ = compare_treatment_contexts(config_a.treatment_context, config_b.treatment_context)
        treat_differs = not treat_match
        mod_match, _ = compare_genetic_modifications(config_a.genetic_modification_context, config_b.genetic_modification_context)
        mod_differs = not mod_match

        if bio_differs or treat_differs or mod_differs:
            differing_dimensions.append(DimensionType.BIOLOGICAL_CONTEXT)
            causes.append(DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT)
            diff_ctx = []
            if bio_differs:
                diff_ctx.append(f"biosample: {config_a.biosample} vs {config_b.biosample}")
            if treat_differs:
                diff_ctx.append(f"treatment: {config_a.treatment_context} vs {config_b.treatment_context}")
            if mod_differs:
                diff_ctx.append(f"genetic_modification: {config_a.genetic_modification_context} vs {config_b.genetic_modification_context}")
            covarying_dims.append("; ".join(diff_ctx))

        # Dimension E: Experiment Assay Discordance (Independent experiments in identical context)
        exp_differs = (
            config_a.encode_experiment_accession
            and config_b.encode_experiment_accession
            and config_a.encode_experiment_accession != config_b.encode_experiment_accession
        )
        if exp_differs and not bio_differs and not treat_differs and not mod_differs:
            differing_dimensions.append(DimensionType.EXPERIMENT_ASSAY)
            causes.append(DisagreementClassification.EXPERIMENT_DISCORDANT)
            covarying_dims.append(
                f"Independent experiments in identical context: {config_a.encode_experiment_accession} vs {config_b.encode_experiment_accession}"
            )

        # 4. Formulate Scientific Explanations and Caveats
        has_multiple_dims = len(differing_dimensions) > 1

        if not causes:
            causes.append(DisagreementClassification.EXPERIMENT_DISCORDANT)
            smallest_diff = "Unspecified experimental or assay variance"
            explanation = "Evidence disagrees across configurations without isolated coordinate or metadata shift."
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
                        f"Promoter window boundary parameter changes exclude or include the peak interval."
                    )
                elif cause == DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT:
                    explanation_parts.append(
                        f"Biological cell line, treatment perturbation, or genetic tag differs between configurations."
                    )
                elif cause == DisagreementClassification.PEAK_PROCESSING_DEPENDENT:
                    explanation_parts.append(
                        f"Peak calling threshold, IDR tier, or replication processing differs within the same experiment."
                    )
                elif cause == DisagreementClassification.EXPERIMENT_DISCORDANT:
                    explanation_parts.append(
                        f"Independent ChIP-seq experiments in the same biological context yield discordant peak detection."
                    )
            explanation = " ".join(explanation_parts)

        # Add mandatory attribution caveat
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
            has_synthetic_configuration=getattr(config_a, "is_synthetic", False) or getattr(config_b, "is_synthetic", False)
        )

    def evaluate_set(
        self,
        configurations: List[EvidenceConfiguration]
    ) -> Dict[str, Any]:
        """
        Evaluates a set of evidence configurations, performing all pairwise comparisons
        and synthesizing the global consensus and disagreement landscape.
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

        for i in range(len(configurations)):
            for j in range(i + 1, len(configurations)):
                diag = self.evaluate_pair(configurations[i], configurations[j])
                # Only record non-stable diagnoses in detailed list
                if DisagreementClassification.STABLE_WITHIN_EVALUATED_CHOICES not in diag.primary_causes:
                    diagnoses.append(diag)
                    for c in diag.primary_causes:
                        observed_causes_set.add(c)

        # Global synthesis
        observed_causes = [c.value for c in sorted(observed_causes_set, key=lambda x: x.value)]

        if not valid_configs:
            overall_status = "INSUFFICIENT_DATA_FAIL_CLOSED"
            interp = (
                f"All {len(configurations)} evaluated configurations were excluded due to unverified "
                f"or missing primary data. Fail-closed policy prevents assigning positive or negative evidence."
            )
        elif len(pos_configs) > 0 and len(neg_configs) == 0:
            overall_status = "ROBUST_POSITIVE_EVIDENCE"
            interp = (
                f"Unanimous positive evidence across {len(pos_configs)} verified configuration(s) "
                f"spanning {len(unique_experiments)} independent experiment(s)."
            )
        elif len(neg_configs) > 0 and len(pos_configs) == 0:
            overall_status = "ROBUST_NO_PEAK_DETECTED"
            interp = (
                f"Unanimous absence of overlapping peaks across {len(neg_configs)} verified configuration(s) "
                f"spanning {len(unique_experiments)} independent experiment(s)."
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
            interp += f" Note: {len(excluded_configs)} configuration(s) were excluded fail-closed due to unverified data."

        is_concordant = bool(len(diagnoses) == 0 and len(valid_configs) > 0)

        return {
            "overall_status": overall_status,
            "total_configurations_count": len(configurations),
            "positive_configurations_count": len(pos_configs),
            "negative_configurations_count": len(neg_configs),
            "excluded_configurations_count": len(excluded_configs),
            "unique_experiments_count": len(unique_experiments),
            "status_concordant": is_concordant,
            "observed_causes": observed_causes,
            "diagnoses": [d.to_dict() for d in diagnoses],
            "scientific_interpretation": interp
        }
