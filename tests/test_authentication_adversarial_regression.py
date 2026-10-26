"""
Adversarial Regression Tests for Source Authentication and Synthetic Separation.
Verifies that manual configurations asserting is_authenticated=True or
auth_status="AUTHENTICATED_PRIMARY_BED" without a sealed ProvenanceToken
are strictly isolated as synthetic mocks and can NEVER achieve
is_verified_biological_observation() == True or contribute to biological consensus.
"""

import pytest
from pathlib import Path
from tf_explorer.disagreement_taxonomy import (
    EvidenceConfiguration,
    EvidenceInterpretation,
    DisagreementClassification,
)
from tf_explorer.disagreement_engine import DisagreementEngine


def test_adversarial_manual_auth_assertion_fails_biological_observation():
    """An adversary creating an EvidenceConfiguration with fake auth flags cannot pass biological verification."""
    fake_config = EvidenceConfiguration(
        config_id="FAKE_CFG_001",
        gene_symbol="TP53",
        transcript_id="ENST00000269305",
        canonical_tss=7676594,
        strand="-",
        window_name="Standard",
        promoter_upstream=2000,
        promoter_downstream=500,
        promoter_start=7676094,
        promoter_end=7678594,
        encode_experiment_accession="ENCSR000FAKE",
        encode_file_accession="ENCFF000FAKE",
        biosample="K562",
        treatment_context="Untreated",
        is_authenticated=True,
        auth_status="AUTHENTICATED_PRIMARY_BED",
        has_overlap=True,
        overlap_bp=200,
        overlap_percentage=100.0,
        evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION,
    )

    # 1. Provenance check must be False
    assert fake_config.has_valid_provenance() is False

    # 2. Synthetic flag must be True (isolated as synthetic mock)
    assert fake_config.is_synthetic is True

    # 3. is_authenticated must be False (cannot claim production authentication)
    assert fake_config.is_authenticated is False

    # 4. auth_status must be SYNTHETIC_DEMONSTRATION, NEVER AUTHENTICATED_PRIMARY_BED
    assert fake_config.auth_status == "SYNTHETIC_DEMONSTRATION"

    # 5. Biological observation verification must be strictly False
    assert fake_config.is_verified_biological_observation() is False

    # 6. Serialized dictionary must reflect synthetic unauthenticated state
    d = fake_config.to_dict()
    assert d["is_synthetic"] is True
    assert d["is_authenticated"] is False
    assert d["auth_status"] == "SYNTHETIC_DEMONSTRATION"


def test_adversarial_synthetic_pair_flags_diagnosis_as_synthetic():
    """DisagreementEngine evaluates geometric properties but strictly marks diagnosis as synthetic."""
    cfg_a = EvidenceConfiguration(
        config_id="SYN_A",
        gene_symbol="TEST",
        transcript_id="ENST0001",
        canonical_tss=10000,
        strand="+",
        window_name="Narrow",
        promoter_upstream=1000,
        promoter_downstream=200,
        promoter_start=9000,
        promoter_end=10200,
        encode_experiment_accession="ENCSR000AAA",
        encode_file_accession="ENCFF000AAA",
        biosample="K562",
        treatment_context="Untreated",
        is_authenticated=False,
        auth_status="SYNTHETIC_DEMONSTRATION",
        is_synthetic=True,
        has_overlap=False,
        evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE,
    )
    cfg_b = EvidenceConfiguration(
        config_id="SYN_B",
        gene_symbol="TEST",
        transcript_id="ENST0001",
        canonical_tss=10000,
        strand="+",
        window_name="Standard",
        promoter_upstream=2000,
        promoter_downstream=500,
        promoter_start=8000,
        promoter_end=10500,
        encode_experiment_accession="ENCSR000AAA",
        encode_file_accession="ENCFF000AAA",
        biosample="K562",
        treatment_context="Untreated",
        is_authenticated=False,
        auth_status="SYNTHETIC_DEMONSTRATION",
        is_synthetic=True,
        has_overlap=True,
        overlap_bp=100,
        evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE,
    )

    engine = DisagreementEngine()
    diag = engine.evaluate_pair(cfg_a, cfg_b, allow_synthetic_fixtures=True)

    # Diagnoses window discordance geometrically in fixture mode
    assert DisagreementClassification.PROMOTER_WINDOW_DEPENDENT in diag.primary_causes

    # Strictly tags diagnosis with synthetic flag
    assert diag.has_synthetic_configuration is True


def test_adversarial_synthetic_concordance_cannot_establish_biological_consensus():
    """Two synthetic configurations claiming identical positive evidence are rejected from consensus."""
    cfg_a = EvidenceConfiguration(
        config_id="SYN_POS_A",
        gene_symbol="TEST",
        transcript_id="ENST0001",
        canonical_tss=10000,
        strand="+",
        window_name="Standard",
        promoter_upstream=2000,
        promoter_downstream=500,
        promoter_start=8000,
        promoter_end=10500,
        encode_experiment_accession="ENCSR000AAA",
        encode_file_accession="ENCFF000AAA",
        biosample="K562",
        treatment_context="Untreated",
        is_authenticated=False,
        auth_status="SYNTHETIC_DEMONSTRATION",
        is_synthetic=True,
        has_overlap=True,
        overlap_bp=100,
        evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE,
    )
    cfg_b = EvidenceConfiguration(
        config_id="SYN_POS_B",
        gene_symbol="TEST",
        transcript_id="ENST0001",
        canonical_tss=10000,
        strand="+",
        window_name="Standard",
        promoter_upstream=2000,
        promoter_downstream=500,
        promoter_start=8000,
        promoter_end=10500,
        encode_experiment_accession="ENCSR000AAA",
        encode_file_accession="ENCFF000AAA",
        biosample="K562",
        treatment_context="Untreated",
        is_authenticated=False,
        auth_status="SYNTHETIC_DEMONSTRATION",
        is_synthetic=True,
        has_overlap=True,
        overlap_bp=100,
        evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE,
    )

    engine = DisagreementEngine()
    diag = engine.evaluate_pair(cfg_a, cfg_b)

    # Cannot establish biological consensus without authentic provenance token
    assert DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA in diag.primary_causes
    assert diag.has_synthetic_configuration is True


def test_synthetic_vs_biological_observation_contract_distinction():
    """
    Verifies that is_verified_observation() is valid within synthetic domain,
    while is_verified_biological_observation() is strictly False.
    """
    cfg = EvidenceConfiguration(
        config_id="SYNTHETIC_DOMAIN_TEST",
        gene_symbol="BRCA1",
        transcript_id="ENST00000357654",
        canonical_tss=43125271,
        strand="-",
        window_name="Standard",
        promoter_upstream=2000,
        promoter_downstream=500,
        promoter_start=43124771,
        promoter_end=43127271,
        is_synthetic=True,
        auth_status="SYNTHETIC_DEMONSTRATION",
        has_overlap=True,
        overlap_bp=150,
        evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE,
    )
    # Valid mathematical synthetic observation:
    assert cfg.is_verified_observation() is True
    # Strictly invalid biological observation:
    assert cfg.is_verified_biological_observation() is False
    assert cfg.is_authenticated is False
    assert cfg.auth_status == "SYNTHETIC_DEMONSTRATION"
