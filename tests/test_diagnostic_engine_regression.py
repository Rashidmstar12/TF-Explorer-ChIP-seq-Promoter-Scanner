"""
Rigorous Regression Test Suite for Patched TF-Explorer Diagnostic Engine.
Validates Phases 2, 3, 4, 5:
- REG-01: Synthetic positive vs synthetic negative (identical params, NEG-E)
- REG-02: Synthetic positive vs synthetic negative (differing window, fixture mode)
- REG-03: Synthetic vs authenticated observation
- REG-04: Unverified vs authenticated observation
- REG-05: Invalid provenance binding (token replay)
- REG-06: Source file mutation (tampered bytes)
- REG-07: Valid authenticated comparison (concordant positive, identical)
- REG-08: Valid authenticated comparison (concordant positive, multi-peak shift)
- REG-09: Valid authenticated comparison (concordant positive, window shift)
- REG-10: Valid authenticated comparison (concordant positive, biosample shift)
- REG-11: Valid authenticated comparison (replicate experiment discordance)
- REG-12: Valid authenticated comparison (concordant negative, no peaks)
"""

import csv
import os
import shutil
from pathlib import Path
import pytest

from tf_explorer.disagreement_taxonomy import (
    DisagreementClassification,
    EvidenceInterpretation,
    DimensionType,
    EvidenceConfiguration,
    DisagreementDiagnosis
)
from tf_explorer.disagreement_engine import DisagreementEngine
from tf_explorer.disagreement_provenance import (
    create_synthetic_provenance_token_for_testing
)
from tf_explorer.disagreement_factory import (
    PrimaryPeakAuthenticator,
    build_validated_configuration
)


CACHE_DIR = (Path(__file__).resolve().parent.parent / "encode_cache").resolve()


@pytest.fixture
def engine():
    return DisagreementEngine()


@pytest.fixture
def auth():
    return PrimaryPeakAuthenticator()


def test_reg01_synthetic_positive_vs_synthetic_negative_identical_params(engine):
    """REG-01: Synthetic positive vs synthetic negative with identical params (NEG-E) must fail closed as INSUFFICIENT_COMPARABLE_DATA, not EXPERIMENT_DISCORDANT."""
    tok_a = create_synthetic_provenance_token_for_testing(file_accession="ENCFF_E1", experiment_accession="ENCSR_SAME", target_tf="E2F1", biosample="MCF-7")
    tok_b = create_synthetic_provenance_token_for_testing(file_accession="ENCFF_E1", experiment_accession="ENCSR_SAME", target_tf="E2F1", biosample="MCF-7")
    cfg_a = EvidenceConfiguration(
        config_id="CFG_E1", gene_symbol="ACTB", canonical_tss=5000, strand="+",
        tf_target="E2F1", encode_experiment_accession="ENCSR_SAME", encode_file_accession="ENCFF_E1",
        biosample="MCF-7", is_authenticated=True, auth_status="AUTHENTICATED_PRIMARY_BED",
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION, provenance_token=tok_a, is_synthetic=True
    )
    cfg_b = EvidenceConfiguration(
        config_id="CFG_E2", gene_symbol="ACTB", canonical_tss=5000, strand="+",
        tf_target="E2F1", encode_experiment_accession="ENCSR_SAME", encode_file_accession="ENCFF_E1",
        biosample="MCF-7", is_authenticated=True, auth_status="AUTHENTICATED_PRIMARY_BED",
        evidence_status=EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED, provenance_token=tok_b, is_synthetic=True
    )
    diag = engine.evaluate_pair(cfg_a, cfg_b)
    assert DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA in diag.primary_causes
    assert DisagreementClassification.EXPERIMENT_DISCORDANT not in diag.primary_causes
    assert "Contradictory evidence statuses asserted on identical experimental accession" in diag.smallest_observed_difference


def test_reg02_synthetic_positive_vs_synthetic_negative_differing_window(engine):
    """REG-02: Synthetic positive vs synthetic negative with window difference."""
    tok_a = create_synthetic_provenance_token_for_testing(file_accession="ENCFF_W1", experiment_accession="ENCSR_W", target_tf="E2F1", biosample="MCF-7")
    tok_b = create_synthetic_provenance_token_for_testing(file_accession="ENCFF_W1", experiment_accession="ENCSR_W", target_tf="E2F1", biosample="MCF-7")
    cfg_a = EvidenceConfiguration(
        config_id="CFG_W1", gene_symbol="ACTB", canonical_tss=5000, strand="+",
        window_start_offset=-2000, window_end_offset=500, promoter_upstream=2000, promoter_downstream=500,
        tf_target="E2F1", encode_experiment_accession="ENCSR_W", encode_file_accession="ENCFF_W1",
        biosample="MCF-7", is_authenticated=True, auth_status="AUTHENTICATED_PRIMARY_BED",
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION, provenance_token=tok_a, is_synthetic=True
    )
    cfg_b = EvidenceConfiguration(
        config_id="CFG_W2", gene_symbol="ACTB", canonical_tss=5000, strand="+",
        window_start_offset=-500, window_end_offset=100, promoter_upstream=500, promoter_downstream=100,
        tf_target="E2F1", encode_experiment_accession="ENCSR_W", encode_file_accession="ENCFF_W1",
        biosample="MCF-7", is_authenticated=True, auth_status="AUTHENTICATED_PRIMARY_BED",
        evidence_status=EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED, provenance_token=tok_b, is_synthetic=True
    )
    diag = engine.evaluate_pair(cfg_a, cfg_b)
    assert DisagreementClassification.PROMOTER_WINDOW_DEPENDENT in diag.primary_causes
    assert diag.has_synthetic_configuration is True


def test_reg03_synthetic_vs_authenticated_observation(engine, auth):
    """REG-03: Comparing synthetic mock against authentic observation must fail closed as INSUFFICIENT_COMPARABLE_DATA."""
    b_encff115 = (CACHE_DIR / "ENCFF115KRO.bed.gz").resolve()
    cfg_auth = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=b_encff115
    )
    tok_syn = create_synthetic_provenance_token_for_testing(file_accession="ENCFF_SYN", experiment_accession="ENCSR_SYN")
    cfg_syn = EvidenceConfiguration(
        config_id="CFG_SYN", gene_symbol="CDK4", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR_SYN", encode_file_accession="ENCFF_SYN",
        biosample="MCF-7", is_authenticated=True, auth_status="AUTHENTICATED_PRIMARY_BED",
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION, provenance_token=tok_syn, is_synthetic=True
    )
    diag = engine.evaluate_pair(cfg_auth, cfg_syn)
    assert DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA in diag.primary_causes
    assert diag.has_synthetic_configuration is True
    assert "authenticated observations cannot be compared against unverified or synthetic" in diag.scientific_explanation


def test_reg04_unverified_vs_authenticated_observation(engine, auth):
    """REG-04: Comparing unverified observation against authenticated observation must fail closed."""
    b_encff115 = (CACHE_DIR / "ENCFF115KRO.bed.gz").resolve()
    cfg_auth = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=b_encff115
    )
    cfg_unauth = EvidenceConfiguration(
        config_id="CFG_UNAUTH", gene_symbol="CDK4", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", is_authenticated=False, auth_status="SOURCE_NOT_VERIFIED",
        evidence_status=EvidenceInterpretation.SOURCE_NOT_VERIFIED
    )
    diag = engine.evaluate_pair(cfg_auth, cfg_unauth)
    assert DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA in diag.primary_causes
    assert "lack authenticated, unmutated primary data" in diag.scientific_explanation


def test_reg05_invalid_provenance_token_binding(engine, auth):
    """REG-05: Replaying a CDK4 token on a MYC configuration must fail closed."""
    b_encff115 = (CACHE_DIR / "ENCFF115KRO.bed.gz").resolve()
    cfg_cdk4 = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=b_encff115
    )
    # Replay attack: attach CDK4 token to MYC configuration
    cfg_replay = EvidenceConfiguration(
        config_id="CFG_REPLAY", gene_symbol="MYC", canonical_tss=127736231, strand="+",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", is_authenticated=True, auth_status="AUTHENTICATED_PRIMARY_BED",
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
        provenance_token=cfg_cdk4.provenance_token, source_file_sha256=cfg_cdk4.source_file_sha256
    )
    diag = engine.evaluate_pair(cfg_cdk4, cfg_replay)
    assert DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA in diag.primary_causes
    assert "Token-configuration binding failed" in diag.smallest_observed_difference


def test_reg06_source_file_mutation(engine, auth, tmp_path):
    """REG-06: Altering bytes in the BED file must fail closed on SHA-256 mismatch."""
    b_encff115 = (CACHE_DIR / "ENCFF115KRO.bed.gz").resolve()
    corrupt_bed = tmp_path / "ENCFF115KRO_mutated.bed.gz"
    shutil.copyfile(b_encff115, corrupt_bed)
    with open(corrupt_bed, "r+b") as f:
        f.seek(50)
        orig = f.read(1)
        f.seek(50)
        f.write(bytes([orig[0] ^ 0xFF]))

    cfg_tampered = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=corrupt_bed
    )
    diag = engine.evaluate_pair(cfg_tampered, cfg_tampered)
    assert DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA in diag.primary_causes
    assert ("Cryptographic hash mismatch" in diag.smallest_observed_difference or
            "Failed to parse source file bytes" in diag.smallest_observed_difference)


def test_reg07_valid_authenticated_concordant_positive_identical(engine, auth):
    """REG-07: Valid authenticated comparison with identical configuration yields concordant positive with configuration concordance."""
    b_encff115 = (CACHE_DIR / "ENCFF115KRO.bed.gz").resolve()
    cfg_a = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=b_encff115
    )
    diag = engine.evaluate_pair(cfg_a, cfg_a)
    assert DisagreementClassification.STABLE_WITHIN_EVALUATED_CHOICES in diag.primary_causes
    assert diag.detection_status_concordance == "BOTH_POSITIVE"
    assert diag.configuration_concordance is True
    assert diag.peak_profile_concordance in ("IDENTICAL_PEAKS", "CONCORDANT_DETECTION")


def test_reg08_valid_authenticated_concordant_positive_multi_peak_shift(engine, auth):
    """REG-08: MYC Standard (2 peaks) vs Distal (1 peak) correctly preserves detection concordance while diagnosing PROMOTER_WINDOW_DEPENDENT and differing peak profiles!"""
    b_encff185 = (CACHE_DIR / "ENCFF185UUT.bed.gz").resolve()
    cfg_std = build_validated_configuration(
        gene_symbol="MYC", transcript_id="ENST00000621592.8", canonical_tss=127736231, strand="+",
        tf_target="ESR1", encode_experiment_accession="ENCSR463GOT", encode_file_accession="ENCFF185UUT",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr8", authenticator=auth, custom_bed_path=b_encff185
    )
    cfg_distal = build_validated_configuration(
        gene_symbol="MYC", transcript_id="ENST00000621592.8", canonical_tss=127736231, strand="+",
        tf_target="ESR1", encode_experiment_accession="ENCSR463GOT", encode_file_accession="ENCFF185UUT",
        biosample="MCF-7", window_name="Distal [-3000, -1000]", window_start_offset=-3000, window_end_offset=-1000,
        chromosome="chr8", authenticator=auth, custom_bed_path=b_encff185
    )
    diag = engine.evaluate_pair(cfg_std, cfg_distal)
    assert diag.detection_status_concordance == "BOTH_POSITIVE"
    assert DisagreementClassification.PROMOTER_WINDOW_DEPENDENT in diag.primary_causes
    assert diag.configuration_concordance is False
    assert "Window definition: [-2000, 500] vs [-3000, -1000]" in diag.smallest_observed_difference


def test_reg09_valid_authenticated_concordant_positive_window_shift(engine, auth):
    """REG-09: CDK4 Standard [-2000, 500] vs Shifted [-1000, 200] (NEG-B) must report PROMOTER_WINDOW_DEPENDENT."""
    b_encff115 = (CACHE_DIR / "ENCFF115KRO.bed.gz").resolve()
    cfg_a = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=b_encff115
    )
    cfg_b = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Shifted [-1000, 200]", window_start_offset=-1000, window_end_offset=200,
        chromosome="chr12", authenticator=auth, custom_bed_path=b_encff115
    )
    diag = engine.evaluate_pair(cfg_a, cfg_b)
    assert diag.detection_status_concordance == "BOTH_POSITIVE"
    assert DisagreementClassification.PROMOTER_WINDOW_DEPENDENT in diag.primary_causes
    assert diag.configuration_concordance is False


def test_reg10_valid_authenticated_concordant_positive_biosample_shift(engine, auth):
    """REG-10: CDK4 E2F1 in MCF-7 vs HeLa-S3 (NEG-C) must report BIOLOGICAL_CONTEXT_DIFFERENT."""
    b_encff115 = (CACHE_DIR / "ENCFF115KRO.bed.gz").resolve()
    b_encff004 = (CACHE_DIR / "ENCFF004ZKN.bed.gz").resolve()
    cfg_a = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EWX", encode_file_accession="ENCFF115KRO",
        biosample="MCF-7", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=b_encff115
    )
    cfg_b = build_validated_configuration(
        gene_symbol="CDK4", transcript_id="ENST00000257904.11", canonical_tss=57752310, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EVM", encode_file_accession="ENCFF004ZKN",
        biosample="HeLa-S3", window_name="Standard [-2000, 500]", window_start_offset=-2000, window_end_offset=500,
        chromosome="chr12", authenticator=auth, custom_bed_path=b_encff004
    )
    diag = engine.evaluate_pair(cfg_a, cfg_b)
    assert diag.detection_status_concordance == "BOTH_POSITIVE"
    assert DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT in diag.primary_causes
    assert diag.configuration_concordance is False


def test_reg11_valid_authenticated_replicate_experiment_discordance(engine, auth):
    """REG-11: PKM Distal CREB1 replicates (SCEN-14) must report EXPERIMENT_DISCORDANT."""
    b_encff142 = (CACHE_DIR / "ENCFF142BHB.bed.gz").resolve()
    b_encff030 = (CACHE_DIR / "ENCFF030DLI.bed.gz").resolve()
    cfg_a = build_validated_configuration(
        gene_symbol="PKM", transcript_id="ENST00000335181.10", canonical_tss=72231190, strand="-",
        tf_target="CREB1", encode_experiment_accession="ENCSR620DUQ", encode_file_accession="ENCFF142BHB",
        biosample="MCF-7", window_name="Distal [-3000, -1000]", window_start_offset=-3000, window_end_offset=-1000,
        chromosome="chr15", authenticator=auth, custom_bed_path=b_encff142
    )
    cfg_b = build_validated_configuration(
        gene_symbol="PKM", transcript_id="ENST00000335181.10", canonical_tss=72231190, strand="-",
        tf_target="CREB1", encode_experiment_accession="ENCSR897JAS", encode_file_accession="ENCFF030DLI",
        biosample="MCF-7", window_name="Distal [-3000, -1000]", window_start_offset=-3000, window_end_offset=-1000,
        chromosome="chr15", authenticator=auth, custom_bed_path=b_encff030
    )
    diag = engine.evaluate_pair(cfg_a, cfg_b)
    assert DisagreementClassification.EXPERIMENT_DISCORDANT in diag.primary_causes
    assert diag.detection_status_concordance == "DISCORDANT"
    assert "Independent ChIP-seq experiments in the same biological context" in diag.scientific_explanation


def test_reg12_valid_authenticated_concordant_negative(engine, auth):
    """REG-12: TFF1 with E2F1 in HeLa-S3 Distal vs Core (both negative) yields concordant non-detection with appropriate disclaimer."""
    b_encff004 = (CACHE_DIR / "ENCFF004ZKN.bed.gz").resolve()
    cfg_a = build_validated_configuration(
        gene_symbol="TFF1", transcript_id="ENST00000291527.3", canonical_tss=42366535, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EVM", encode_file_accession="ENCFF004ZKN",
        biosample="HeLa-S3", window_name="Distal [-3000, -1000]", window_start_offset=-3000, window_end_offset=-1000,
        chromosome="chr21", authenticator=auth, custom_bed_path=b_encff004
    )
    cfg_b = build_validated_configuration(
        gene_symbol="TFF1", transcript_id="ENST00000291527.3", canonical_tss=42366535, strand="-",
        tf_target="E2F1", encode_experiment_accession="ENCSR000EVM", encode_file_accession="ENCFF004ZKN",
        biosample="HeLa-S3", window_name="Core [-1000, 100]", window_start_offset=-1000, window_end_offset=100,
        chromosome="chr21", authenticator=auth, custom_bed_path=b_encff004
    )
    diag = engine.evaluate_pair(cfg_a, cfg_b)
    assert diag.detection_status_concordance == "BOTH_NEGATIVE"
    assert diag.peak_profile_concordance == "NO_PEAKS"
    assert "does not establish biological non-binding across unassayed contexts" in diag.causal_attribution_caveat
