"""
Quickstart Example: Regulatory Disagreement Engine
Demonstrates automated discordance triage and classification of conflicting regulatory evidence.
"""

from tf_explorer.disagreement_taxonomy import (
    DisagreementClassification,
    EvidenceInterpretation,
    EvidenceConfiguration,
)
from tf_explorer.disagreement_engine import DisagreementEngine

def main():
    print("=== TF-Explorer: Disagreement Engine Demonstration ===")
    
    engine = DisagreementEngine()
    
    # Scenario: Peak observed at -1500 bp from TSS of TP53.
    # Narrow window [-1000, +200] excludes peak -> NEGATIVE_EVIDENCE
    # Standard window [-2000, +500] includes peak -> POSITIVE_EVIDENCE
    cfg_narrow = EvidenceConfiguration(
        config_id="CFG_NARROW_WINDOW",
        gene_symbol="TP53",
        transcript_id="ENST00000269305",
        canonical_tss=7676594,
        strand="-",
        window_name="Narrow",
        promoter_upstream=1000,
        promoter_downstream=200,
        promoter_start=7676394,
        promoter_end=7677594,
        encode_experiment_accession="ENCSR000BKS",
        encode_file_accession="ENCFF001AAA",
        biosample="K562",
        treatment_context="Untreated",
        is_authenticated=False,
        auth_status="SYNTHETIC_DEMONSTRATION",
        is_synthetic=True,
        synthetic_note="Synthetic educational demonstration of promoter window discordance",
        has_overlap=False,
        overlap_bp=0,
        overlap_percentage=0.0,
        evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE,
    )
    
    cfg_standard = EvidenceConfiguration(
        config_id="CFG_STANDARD_WINDOW",
        gene_symbol="TP53",
        transcript_id="ENST00000269305",
        canonical_tss=7676594,
        strand="-",
        window_name="Standard",
        promoter_upstream=2000,
        promoter_downstream=500,
        promoter_start=7676094,
        promoter_end=7678594,
        encode_experiment_accession="ENCSR000BKS",
        encode_file_accession="ENCFF001AAA",
        biosample="K562",
        treatment_context="Untreated",
        is_authenticated=False,
        auth_status="SYNTHETIC_DEMONSTRATION",
        is_synthetic=True,
        synthetic_note="Synthetic educational demonstration of promoter window discordance",
        has_overlap=True,
        overlap_bp=150,
        overlap_percentage=100.0,
        evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE,
    )
    
    diagnosis = engine.evaluate_pair(cfg_narrow, cfg_standard)
    
    print(f"Narrow promoter window [-1000, +200]: {cfg_narrow.evidence_status.value}")
    print(f"Standard promoter window [-2000, +500]: {cfg_standard.evidence_status.value}")
    print(f"Automated Dispute Taxonomy Classification: {[c.value for c in diagnosis.primary_causes]}")
    print(f"Differing Dimensions: {[d.value for d in diagnosis.differing_dimensions]}")
    print(f"Scientific Explanation: {diagnosis.scientific_explanation}")
    print(f"Synthetic Demonstration Flag: {diagnosis.has_synthetic_configuration}")
    
    assert DisagreementClassification.PROMOTER_WINDOW_DEPENDENT in diagnosis.primary_causes
    assert diagnosis.has_synthetic_configuration is True
    assert cfg_narrow.is_verified_biological_observation() is False
    assert cfg_standard.is_verified_biological_observation() is False
    assert cfg_narrow.is_authenticated is False
    assert cfg_standard.is_authenticated is False
    assert cfg_narrow.auth_status == "SYNTHETIC_DEMONSTRATION"
    assert cfg_standard.auth_status == "SYNTHETIC_DEMONSTRATION"
    print("\n[PROVENANCE GUARD] Synthetic configuration verified: strictly unauthenticated and excluded from biological consensus.")
    print("[SUCCESS] Disagreement Engine successfully classified regulatory discordance in synthetic demonstration.")

if __name__ == "__main__":
    main()
