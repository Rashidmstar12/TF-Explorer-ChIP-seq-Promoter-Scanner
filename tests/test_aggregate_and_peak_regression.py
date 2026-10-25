"""
Tests for evaluate_set() aggregate analysis and peak-profile comparability.
Covers Cases A through H and full physical peak comparability rules.
"""

import pytest
from pathlib import Path

from tf_explorer.disagreement_taxonomy import (
    EvidenceConfiguration,
    EvidenceInterpretation,
    DisagreementClassification,
    DimensionType
)
from tf_explorer.disagreement_engine import (
    DisagreementEngine,
    classify_peak_profile_concordance
)
from tf_explorer.disagreement_factory import (
    PrimaryPeakAuthenticator,
    build_validated_configuration
)
from tf_explorer.disagreement_provenance import (
    create_synthetic_provenance_token_for_testing
)

CACHE_DIR = (Path(__file__).resolve().parent.parent / "encode_cache").resolve()
B_ENCFF115 = CACHE_DIR / "ENCFF115KRO.bed.gz"
B_ENCFF004 = CACHE_DIR / "ENCFF004ZKN.bed.gz"
B_ENCFF185 = CACHE_DIR / "ENCFF185UUT.bed.gz"


@pytest.fixture
def engine():
    return DisagreementEngine()


@pytest.fixture
def auth():
    return PrimaryPeakAuthenticator()


# =========================================================================
# Section 2: Cases A through H (Aggregate Analysis Tests)
# =========================================================================

def test_case_a_differing_promoter_windows(engine, auth):
    """
    Case A: Two positive configurations with different promoter windows.
    Verifies that evaluate_set() retains the PROMOTER_WINDOW_DEPENDENT diagnosis
    and does not silently discard it.
    """
    cfg_a = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=B_ENCFF115
    )
    cfg_b = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Shifted [-1000, 200]", window_start_offset=-1000, window_end_offset=200,
        chromosome="chr12", authenticator=auth, custom_bed_path=B_ENCFF115
    )

    # Pairwise check
    diag = engine.evaluate_pair(cfg_a, cfg_b)
    assert diag.detection_status_concordance == "BOTH_POSITIVE"
    assert diag.configuration_concordance is False
    assert DisagreementClassification.PROMOTER_WINDOW_DEPENDENT in diag.primary_causes

    # Aggregate check
    res = engine.evaluate_set([cfg_a, cfg_b])
    assert res["positive_configurations_count"] == 2
    assert res["negative_configurations_count"] == 0
    assert res["status_concordant"] is True
    assert res["detection_status_concordant"] is True
    assert res["configuration_concordant"] is False
    assert "PROMOTER_WINDOW_DEPENDENT" in res["observed_causes"]
    assert len(res["diagnoses"]) == 1
    assert res["overall_status"] == "CONCORDANT_POSITIVE_PARAMETER_SENSITIVE"
    assert "PROMOTER_WINDOW_DEPENDENT" in res["scientific_interpretation"]
    assert "Observations must not be assumed to have identical genomic boundaries" in res["scientific_interpretation"]


def test_case_b_differing_physical_peaks(engine, auth):
    """
    Case B: Two positive configurations with different physical peak identities.
    (MYC ESR1 in MCF-7: Standard window captures 2 peaks, Distal captures 1 peak).
    """
    cfg_a = build_validated_configuration(
        gene_symbol="MYC", transcript_id="ENST00000621592.8", canonical_tss=127736231, strand="+",
        tf_target="ESR1", encode_experiment_accession="ENCSR463GOT", encode_file_accession="ENCFF185UUT",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr8", authenticator=auth, custom_bed_path=B_ENCFF185
    )
    cfg_b = build_validated_configuration(
        gene_symbol="MYC", transcript_id="ENST00000621592.8", canonical_tss=127736231, strand="+",
        tf_target="ESR1", encode_experiment_accession="ENCSR463GOT", encode_file_accession="ENCFF185UUT",
        biosample="MCF-7", window_name="Distal [-3000, -1000]", window_start_offset=-3000, window_end_offset=-1000,
        chromosome="chr8", authenticator=auth, custom_bed_path=B_ENCFF185
    )

    # Pairwise check
    diag = engine.evaluate_pair(cfg_a, cfg_b)
    assert diag.detection_status_concordance == "BOTH_POSITIVE"
    assert diag.peak_profile_concordance == "DIFFERING_PEAK_COUNTS"

    # Aggregate check
    res = engine.evaluate_set([cfg_a, cfg_b])
    assert res["positive_configurations_count"] == 2
    assert res["status_concordant"] is True
    assert res["peak_profile_concordant"] is False
    assert len(res["diagnoses"]) == 1
    assert "DIFFERING_PEAK_COUNTS" in res["diagnoses"][0]["peak_profile_concordance"]
    assert "Observations must not be assumed to have identical genomic boundaries or peak multiplicity" in res["scientific_interpretation"]


def test_case_c_differing_biosamples(engine, auth):
    """
    Case C: Two positive configurations from different biosamples (MCF-7 vs HeLa-S3).
    """
    cfg_a = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=B_ENCFF115
    )
    cfg_b = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EVM", encode_file_accession="ENCFF004ZKN",
        biosample="HeLa-S3", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=B_ENCFF004
    )

    # Pairwise check
    diag = engine.evaluate_pair(cfg_a, cfg_b)
    assert diag.detection_status_concordance == "BOTH_POSITIVE"
    assert DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT in diag.primary_causes

    # Aggregate check
    res = engine.evaluate_set([cfg_a, cfg_b])
    assert res["positive_configurations_count"] == 2
    assert res["status_concordant"] is True
    assert res["configuration_concordant"] is False
    assert "BIOLOGICAL_CONTEXT_DIFFERENT" in res["observed_causes"]
    assert len(res["diagnoses"]) == 1
    assert res["overall_status"] == "CONCORDANT_POSITIVE_PARAMETER_SENSITIVE"


def test_case_d_mixed_positive_and_negative(engine, auth):
    """
    Case D: Three configurations with mixed positive, positive, and negative observations.
    """
    cfg_pos1 = build_validated_configuration(
        gene_symbol="TFF1", transcript_id="ENST00000291527.3", canonical_tss=42366535, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr21", authenticator=auth, custom_bed_path=B_ENCFF115
    )
    cfg_pos2 = build_validated_configuration(
        gene_symbol="TFF1", transcript_id="ENST00000291527.3", canonical_tss=42366535, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Shifted [-1000, 200]", window_start_offset=-1000, window_end_offset=200,
        chromosome="chr21", authenticator=auth, custom_bed_path=B_ENCFF115
    )
    cfg_neg = build_validated_configuration(
        gene_symbol="TFF1", transcript_id="ENST00000291527.3", canonical_tss=42366535, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EVM", encode_file_accession="ENCFF004ZKN",
        biosample="HeLa-S3", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr21", authenticator=auth, custom_bed_path=B_ENCFF004
    )

    res = engine.evaluate_set([cfg_pos1, cfg_pos2, cfg_neg])
    assert res["total_configurations_count"] == 3
    assert res["positive_configurations_count"] == 2
    assert res["negative_configurations_count"] == 1
    assert res["status_concordant"] is False
    assert res["detection_status_concordant"] is False
    assert res["overall_status"] == "DISCORDANT_EVIDENCE"
    assert "BIOLOGICAL_CONTEXT_DIFFERENT" in res["observed_causes"]
    assert "PROMOTER_WINDOW_DEPENDENT" in res["observed_causes"]
    assert len(res["diagnoses"]) == 3  # All 3 pairs have differences


def test_case_e_distinct_negative_contexts(engine, auth):
    """
    Case E: Two negative configurations from distinct experimental contexts.
    Verifies that concordant negative is flagged with explicit caveats against claiming non-inactivity.
    """
    cfg_neg1 = build_validated_configuration(
        gene_symbol="TFF1", transcript_id="ENST00000291527.3", canonical_tss=42366535, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EVM", encode_file_accession="ENCFF004ZKN",
        biosample="HeLa-S3", window_name="Core [-500, 100]", window_start_offset=-500, window_end_offset=100,
        chromosome="chr21", authenticator=auth, custom_bed_path=B_ENCFF004
    )
    cfg_neg2 = build_validated_configuration(
        gene_symbol="TFF1", transcript_id="ENST00000291527.3", canonical_tss=42366535, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EVM", encode_file_accession="ENCFF004ZKN",
        biosample="HeLa-S3", window_name="Distal [-3000, -1000]", window_start_offset=-3000, window_end_offset=-1000,
        chromosome="chr21", authenticator=auth, custom_bed_path=B_ENCFF004
    )

    # Pairwise check
    diag = engine.evaluate_pair(cfg_neg1, cfg_neg2)
    assert diag.detection_status_concordance == "BOTH_NEGATIVE"
    assert diag.peak_profile_concordance == "NO_PEAKS"
    assert "does not establish biological non-binding" in diag.causal_attribution_caveat

    # Aggregate check
    res = engine.evaluate_set([cfg_neg1, cfg_neg2])
    assert res["positive_configurations_count"] == 0
    assert res["negative_configurations_count"] == 2
    assert res["status_concordant"] is True
    assert res["overall_status"] == "CONCORDANT_NEGATIVE_PARAMETER_SENSITIVE"
    assert "Non-detection within evaluated promoter windows does not establish biological inactivity" in res["scientific_interpretation"]


def test_case_f_authenticated_and_synthetic(engine, auth):
    """
    Case F: Authenticated and synthetic configurations in the same analysis set.
    Synthetic configuration must be excluded from valid count fail-closed,
    and pairwise diagnosis must report INSUFFICIENT_COMPARABLE_DATA.
    """
    cfg_auth = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", chromosome="chr12", authenticator=auth, custom_bed_path=B_ENCFF115
    )
    tok_syn = create_synthetic_provenance_token_for_testing(
        file_accession="ENCFF_SYN", experiment_accession="ENCSR_SYN", target_tf="E2F1", biosample="MCF-7"
    )
    cfg_syn = EvidenceConfiguration(
        config_id="CFG_SYN", gene_symbol="CDK4", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR_SYN", encode_file_accession="ENCFF_SYN",
        biosample="MCF-7", is_authenticated=True, is_synthetic=True, auth_status="AUTHENTICATED_PRIMARY_BED",
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION, provenance_token=tok_syn,
    )

    res = engine.evaluate_set([cfg_auth, cfg_syn])
    assert res["total_configurations_count"] == 2
    assert res["positive_configurations_count"] == 1
    assert res["excluded_configurations_count"] == 1
    assert res["has_invalid_or_unauthenticated_evidence"] is True
    assert "INSUFFICIENT_COMPARABLE_DATA" in res["observed_causes"]
    assert len(res["diagnoses"]) == 1
    assert "excluded fail-closed" in res["scientific_interpretation"]


def test_case_g_multiple_simultaneous_differences(engine, auth):
    """
    Case G: Multiple simultaneous analytical differences (biosample + window).
    """
    cfg_a = build_validated_configuration(
        gene_symbol="TFF1", transcript_id="ENST00000291527.3", canonical_tss=42366535, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr21", authenticator=auth, custom_bed_path=B_ENCFF115
    )
    cfg_b = build_validated_configuration(
        gene_symbol="TFF1", transcript_id="ENST00000291527.3", canonical_tss=42366535, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EVM", encode_file_accession="ENCFF004ZKN",
        biosample="HeLa-S3", window_name="Distal [-3000, -1000]", window_start_offset=-3000, window_end_offset=-1000,
        chromosome="chr21", authenticator=auth, custom_bed_path=B_ENCFF004
    )

    res = engine.evaluate_set([cfg_a, cfg_b])
    assert "BIOLOGICAL_CONTEXT_DIFFERENT" in res["observed_causes"]
    assert "PROMOTER_WINDOW_DEPENDENT" in res["observed_causes"]
    assert len(res["diagnoses"]) == 1
    diag_dict = res["diagnoses"][0]
    assert diag_dict["has_multiple_covarying_dimensions"] is True
    assert "MULTIPLE_COVARYING_DIMENSIONS" in diag_dict["causal_attribution_caveat"]


def test_case_h_empty_configuration_set(engine):
    """
    Case H: Empty configuration set.
    """
    res = engine.evaluate_set([])
    assert res["overall_status"] == "EMPTY_SET"
    assert res["total_configurations_count"] == 0
    assert res["positive_configurations_count"] == 0
    assert res["negative_configurations_count"] == 0
    assert res["excluded_configurations_count"] == 0
    assert res["status_concordant"] is True
    assert res["observed_causes"] == []
    assert res["diagnoses"] == []
    assert res["scientific_interpretation"] == "No configurations provided for evaluation."


# =========================================================================
# Section 7: Peak-Profile Comparability Unit Tests
# =========================================================================

def test_peak_profile_identical():
    c_a = EvidenceConfiguration(
        config_id="A", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}]
    )
    c_b = EvidenceConfiguration(
        config_id="B", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}]
    )
    res = classify_peak_profile_concordance(c_a, c_b, True, False)
    assert res == "IDENTICAL_PEAKS"


def test_peak_profile_partially_shared():
    c_a = EvidenceConfiguration(
        config_id="A", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}, {"chrom": "chr1", "start": 300, "end": 400}]
    )
    c_b = EvidenceConfiguration(
        config_id="B", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}, {"chrom": "chr1", "start": 500, "end": 600}]
    )
    res = classify_peak_profile_concordance(c_a, c_b, True, False)
    assert res == "PARTIALLY_SHARED_PEAKS"


def test_peak_profile_disjoint():
    c_a = EvidenceConfiguration(
        config_id="A", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}]
    )
    c_b = EvidenceConfiguration(
        config_id="B", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 500, "end": 600}]
    )
    res = classify_peak_profile_concordance(c_a, c_b, True, False)
    assert res == "DISTINCT_PEAKS"


def test_peak_profile_differing_counts():
    c_a = EvidenceConfiguration(
        config_id="A", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}, {"chrom": "chr1", "start": 300, "end": 400}]
    )
    c_b = EvidenceConfiguration(
        config_id="B", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}]
    )
    res = classify_peak_profile_concordance(c_a, c_b, True, False)
    assert res == "DIFFERING_PEAK_COUNTS"


def test_peak_profile_overlapping_offset():
    c_a = EvidenceConfiguration(
        config_id="A", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 250}]
    )
    c_b = EvidenceConfiguration(
        config_id="B", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 120, "end": 270}]
    )
    res = classify_peak_profile_concordance(c_a, c_b, True, False)
    assert res == "OVERLAPPING_PEAKS_OFFSET"


def test_peak_profile_missing_details():
    c_a = EvidenceConfiguration(
        config_id="A", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[], source_peak_start=None, source_peak_end=None
    )
    c_b = EvidenceConfiguration(
        config_id="B", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[], source_peak_start=None, source_peak_end=None
    )
    res = classify_peak_profile_concordance(c_a, c_b, True, False)
    assert res == "INCOMPARABLE"


def test_peak_profile_incompatible_assembly():
    c_a = EvidenceConfiguration(
        config_id="A", gene_symbol="GENE", tf_target="E2F1", assembly="hg19", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}]
    )
    c_b = EvidenceConfiguration(
        config_id="B", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}]
    )
    res = classify_peak_profile_concordance(c_a, c_b, True, False)
    assert res == "INCOMPARABLE"


def test_peak_profile_different_target_tfs():
    c_a = EvidenceConfiguration(
        config_id="A", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}]
    )
    c_b = EvidenceConfiguration(
        config_id="B", gene_symbol="GENE", tf_target="MYC", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        all_overlapping_peaks=[{"chrom": "chr1", "start": 100, "end": 200}]
    )
    res = classify_peak_profile_concordance(c_a, c_b, True, False)
    assert res == "INCOMPARABLE"


def test_peak_profile_concordant_negative():
    c_a = EvidenceConfiguration(
        config_id="A", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED
    )
    c_b = EvidenceConfiguration(
        config_id="B", gene_symbol="GENE", tf_target="E2F1", assembly="GRCh38", is_authenticated=True, is_synthetic=True,
        evidence_status=EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED
    )
    res = classify_peak_profile_concordance(c_a, c_b, False, True)
    assert res == "NO_PEAKS"