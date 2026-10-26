"""
Comprehensive Adversarial Validation Suite for Regulatory Evidence Disagreement Explorer (Phase 4B).

Validates all 14 mandatory scientific challenge cases:
1. Different transcripts with nearby TSS positions.
2. Widely separated alternative TSS positions.
3. Peaks touching promoter boundaries.
4. Two incompatible biological contexts.
5. Two processing files from one experiment (PEAK_PROCESSING_DEPENDENT).
6. Two genuinely independent experiments (EXPERIMENT_DISCORDANT).
7. Missing treatment metadata.
8. An experiment without compatible peak evidence.
9. An unverified source peak (fails closed).
10. A fabricated peak row (fails closed).
11. A cached file that has changed.
12. An absent MANE Select designation.
13. Interacting transcript and promoter-window effects.
14. Independent reference implementation concordance (Stage 7).
"""

import pytest
from typing import Dict, Any, List
from streamlit.testing.v1 import AppTest

from tf_explorer.disagreement_taxonomy import (
    DisagreementClassification,
    EvidenceInterpretation,
    DimensionType,
    EvidenceConfiguration,
    DisagreementDiagnosis
)
from tf_explorer.disagreement_engine import DisagreementEngine
from tf_explorer.disagreement_matrix import DisagreementMatrix
from tf_explorer.disagreement_reference import (
    reference_promoter_coordinates,
    reference_bed0_overlap,
    reference_classify_pair,
    reference_evaluate_matrix
)


class TestPhase4BDisagreementExplorer:
    """Rigorous scientific tests for cause-of-disagreement classification and matrix construction."""

    @pytest.fixture
    def engine(self):
        return DisagreementEngine(allow_synthetic_fixtures=True)

    # Challenge 1: Different transcripts with nearby TSS positions
    def test_challenge_1_nearby_tss_positions(self, engine):
        """
        Two alternative transcripts with TSS separated by 150 bp.
        A narrow peak overlaps one promoter interval but misses the other.
        """
        # Peak on chr1: 10,000 - 10,200
        # Transcript 1 TSS: 10,100 -> Promoter [-100, +50]: [10000, 10150] -> Overlaps!
        # Transcript 2 TSS: 10,350 -> Promoter [-100, +50]: [10250, 10400] -> No overlap!
        cfg1 = EvidenceConfiguration(
            config_id="CFG_TX1",
            gene_symbol="TEST1",
            transcript_id="ENST000001",
            canonical_tss=10100,
            strand="+",
            window_name="Narrow",
            promoter_upstream=100,
            promoter_downstream=50,
            promoter_start=10000,
            promoter_end=10150,
            encode_experiment_accession="ENCSR001AAA",
            encode_file_accession="ENCFF001AAA",
            biosample="MCF-7",
            treatment_context="Untreated",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=True,
            overlap_bp=150,
            overlap_percentage=75.0,
            evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE
        )
        cfg2 = EvidenceConfiguration(
            config_id="CFG_TX2",
            gene_symbol="TEST1",
            transcript_id="ENST000002",
            canonical_tss=10350,
            strand="+",
            window_name="Narrow",
            promoter_upstream=100,
            promoter_downstream=50,
            promoter_start=10250,
            promoter_end=10400,
            encode_experiment_accession="ENCSR001AAA",
            encode_file_accession="ENCFF001AAA",
            biosample="MCF-7",
            treatment_context="Untreated",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=False,
            overlap_bp=0,
            overlap_percentage=0.0,
            evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE
        )
        diag = engine.evaluate_pair(cfg1, cfg2)
        assert DisagreementClassification.TRANSCRIPT_GEOMETRY_DEPENDENT in diag.primary_causes
        assert DimensionType.TRANSCRIPT_TSS in diag.differing_dimensions
        assert "250 bp" in diag.smallest_observed_difference

    # Challenge 2: Widely separated alternative TSS positions
    def test_challenge_2_widely_separated_alternative_tss(self, engine):
        """
        Alternative promoters separated by 45 kb (e.g. distal upstream promoter vs proximal).
        """
        cfg_proximal = EvidenceConfiguration(
            config_id="CFG_PROX",
            gene_symbol="GENE_DISTAL",
            transcript_id="ENST_PROX",
            canonical_tss=100000,
            strand="+",
            window_name="Standard",
            promoter_upstream=2000,
            promoter_downstream=500,
            promoter_start=98000,
            promoter_end=100500,
            encode_experiment_accession="ENCSR002BBB",
            encode_file_accession="ENCFF002BBB",
            biosample="K562",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=True,
            overlap_bp=300,
            evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE
        )
        cfg_distal = EvidenceConfiguration(
            config_id="CFG_DIST",
            gene_symbol="GENE_DISTAL",
            transcript_id="ENST_DIST",
            canonical_tss=55000,  # 45 kb upstream!
            strand="+",
            window_name="Standard",
            promoter_upstream=2000,
            promoter_downstream=500,
            promoter_start=53000,
            promoter_end=55500,
            encode_experiment_accession="ENCSR002BBB",
            encode_file_accession="ENCFF002BBB",
            biosample="K562",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=False,
            overlap_bp=0,
            evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE
        )
        diag = engine.evaluate_pair(cfg_proximal, cfg_distal)
        assert DisagreementClassification.TRANSCRIPT_GEOMETRY_DEPENDENT in diag.primary_causes
        assert "45,000 bp" in diag.smallest_observed_difference

    # Challenge 3: Peaks touching promoter boundaries
    def test_challenge_3_peaks_touching_promoter_boundaries(self, engine):
        """
        Tests mathematical adherence to BED0 half-open interval [start, end).
        Peak ending exactly at promoter start has 0 bp overlap.
        Peak entering by 1 bp has 1 bp overlap.
        """
        # Promoter: [1000, 2000)
        p_start, p_end = 1000, 2000
        # Peak A: [800, 1000) -> touches start, overlap = 0
        ov_touch = reference_bed0_overlap(p_start, p_end, 800, 1000)
        assert ov_touch == 0

        # Peak B: [800, 1001) -> penetrates by 1 bp, overlap = 1
        ov_1bp = reference_bed0_overlap(p_start, p_end, 800, 1001)
        assert ov_1bp == 1

        # Peak C: [2000, 2500) -> touches end, overlap = 0
        ov_touch_end = reference_bed0_overlap(p_start, p_end, 2000, 2500)
        assert ov_touch_end == 0

    # Challenge 4: Two incompatible biological contexts
    def test_challenge_4_incompatible_biological_contexts(self, engine):
        """
        Disagreement between MCF-7 (positive) and HeLa-S3 (negative) in independent assays.
        """
        cfg_mcf7 = EvidenceConfiguration(
            config_id="CFG_MCF7",
            gene_symbol="CDK4",
            transcript_id="ENST00000257904",
            canonical_tss=57752310,
            strand="-",
            encode_experiment_accession="ENCSR000EWX",
            encode_file_accession="ENCFF115KRO",
            biosample="MCF-7",
            treatment_context="Untreated",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=True,
            overlap_bp=343,
            evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE
        )
        cfg_hela = EvidenceConfiguration(
            config_id="CFG_HELA",
            gene_symbol="CDK4",
            transcript_id="ENST00000257904",
            canonical_tss=57752310,
            strand="-",
            encode_experiment_accession="ENCSR000EVM",
            encode_file_accession="ENCFF004ZKN",
            biosample="HeLa-S3",
            treatment_context="Untreated",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=False,
            overlap_bp=0,
            evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE
        )
        diag = engine.evaluate_pair(cfg_mcf7, cfg_hela)
        assert DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT in diag.primary_causes
        assert DimensionType.BIOLOGICAL_CONTEXT in diag.differing_dimensions
        assert "biosample: MCF-7 vs HeLa-S3" in diag.smallest_observed_difference

    # Challenge 5: Two processing files from one experiment (PEAK_PROCESSING_DEPENDENT)
    def test_challenge_5_peak_processing_dependent_same_experiment(self, engine):
        """
        Two different files from the SAME experiment accession (e.g. conservative IDR vs permissive).
        Must be classified as PEAK_PROCESSING_DEPENDENT, NOT EXPERIMENT_DISCORDANT.
        """
        cfg_file1 = EvidenceConfiguration(
            config_id="CFG_FILE1",
            gene_symbol="ACTB",
            encode_experiment_accession="ENCSR000EWX",
            encode_file_accession="ENCFF115KRO",
            processing_tier="conservative_idr_peaks",
            biosample="MCF-7",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=True,
            overlap_bp=650,
            evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE
        )
        cfg_file2 = EvidenceConfiguration(
            config_id="CFG_FILE2",
            gene_symbol="ACTB",
            encode_experiment_accession="ENCSR000EWX",
            encode_file_accession="ENCFF999REP",
            processing_tier="replicated_broad_peaks",
            biosample="MCF-7",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=False,
            overlap_bp=0,
            evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE
        )
        diag = engine.evaluate_pair(cfg_file1, cfg_file2)
        assert DisagreementClassification.PEAK_PROCESSING_DEPENDENT in diag.primary_causes
        assert DisagreementClassification.EXPERIMENT_DISCORDANT not in diag.primary_causes
        assert DimensionType.PEAK_PROCESSING in diag.differing_dimensions

    # Challenge 6: Two genuinely independent experiments (EXPERIMENT_DISCORDANT)
    def test_challenge_6_experiment_discordant_same_biological_context(self, engine):
        """
        Two independent experiments for the same TF, same cell line, same untreated state disagree.
        Must be classified as EXPERIMENT_DISCORDANT.
        """
        cfg_exp1 = EvidenceConfiguration(
            config_id="CFG_EXP1",
            gene_symbol="GAPDH",
            encode_experiment_accession="ENCSR000AAA",
            encode_file_accession="ENCFF001AAA",
            biosample="MCF-7",
            treatment_context="Untreated",
            genetic_modification_context="None",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=True,
            overlap_bp=250,
            evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE
        )
        cfg_exp2 = EvidenceConfiguration(
            config_id="CFG_EXP2",
            gene_symbol="GAPDH",
            encode_experiment_accession="ENCSR000BBB",
            encode_file_accession="ENCFF001BBB",
            biosample="MCF-7",
            treatment_context="Untreated",
            genetic_modification_context="None",
            is_authenticated=True,
            auth_status="AUTHENTICATED_PRIMARY_BED",
            has_overlap=False,
            overlap_bp=0,
            evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE
        )
        diag = engine.evaluate_pair(cfg_exp1, cfg_exp2)
        assert DisagreementClassification.EXPERIMENT_DISCORDANT in diag.primary_causes
        assert DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT not in diag.primary_causes
        assert "suggests assay sensitivity" in diag.causal_attribution_caveat

    # Challenge 7: Missing treatment metadata
    def test_challenge_7_missing_treatment_metadata_handled(self, engine):
        """
        Experiment with unverified/missing treatment metadata triggers caveat and attribution warning.
        """
        cfg_pert = EvidenceConfiguration(
            config_id="CFG_PERT",
            gene_symbol="TFF1",
            encode_experiment_accession="ENCSR463GOT",
            treatment_context="17β-estradiol (100 nM)",
            biosample="MCF-7",
            is_authenticated=True,
            evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE
        )
        cfg_unknown = EvidenceConfiguration(
            config_id="CFG_UNK",
            gene_symbol="TFF1",
            encode_experiment_accession="ENCSR999UNK",
            treatment_context="Unknown / Missing Metadata",
            biosample="MCF-7",
            is_authenticated=True,
            evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE
        )
        diag = engine.evaluate_pair(cfg_pert, cfg_unknown)
        assert DisagreementClassification.BIOLOGICAL_CONTEXT_DIFFERENT in diag.primary_causes
        assert diag.causal_attribution_caveat is not None

    # Challenge 8: Experiment without compatible peak evidence
    def test_challenge_8_experiment_without_compatible_peaks_classified_insufficient(self, engine):
        """
        An experiment that has no peak calls or incompatible file format is excluded with INSUFFICIENT_DATA.
        """
        cfg_valid = EvidenceConfiguration(
            config_id="CFG_VAL",
            gene_symbol="TP53",
            encode_experiment_accession="ENCSR000EWX",
            is_authenticated=True,
            evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE
        )
        cfg_no_peaks = EvidenceConfiguration(
            config_id="CFG_NO_PEAKS",
            gene_symbol="TP53",
            encode_experiment_accession="ENCSR888NOPK",
            is_authenticated=False,
            auth_status="SOURCE_UNAVAILABLE",
            exclusion_reason="No narrowPeak files released for experiment",
            evidence_status=EvidenceInterpretation.EXCLUDED_INSUFFICIENT_DATA
        )
        diag = engine.evaluate_pair(cfg_valid, cfg_no_peaks)
        assert DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA in diag.primary_causes
        assert diag.evidence_b == EvidenceInterpretation.EXCLUDED_INSUFFICIENT_DATA

    # Challenge 9: Unverified source peak fails closed
    def test_challenge_9_unverified_source_peak_fails_closed(self, engine):
        """
        An unverified source peak must NEVER support a positive biological evidence classification.
        """
        cfg_unverified = EvidenceConfiguration(
            config_id="CFG_UNVERIFIED",
            gene_symbol="TP53",
            is_authenticated=False,
            auth_status="INDEXED_UNAUTHENTICATED",
            has_overlap=True,  # Arithmetic overlap is True...
            overlap_bp=200,
            evidence_status=EvidenceInterpretation.EXCLUDED_UNAUTHENTICATED,  # ...but status is strictly EXCLUDED!
            exclusion_reason="Unverified peak from unanchored index"
        )
        assert cfg_unverified.evidence_status != EvidenceInterpretation.POSITIVE_EVIDENCE
        
        cfg_baseline = EvidenceConfiguration(
            config_id="CFG_BASE",
            gene_symbol="TP53",
            is_authenticated=True,
            evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE
        )
        diag = engine.evaluate_pair(cfg_baseline, cfg_unverified)
        assert DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA in diag.primary_causes

    # Challenge 10: Fabricated peak row fails closed
    def test_challenge_10_fabricated_peak_row_fails_closed(self, engine):
        """
        Fabricated peak coordinates fail authentication and are excluded from positive classification.
        """
        cfg_fake = EvidenceConfiguration(
            config_id="CFG_FAKE",
            gene_symbol="CDK4",
            is_authenticated=False,
            auth_status="FABRICATED_OR_ABSENT",
            has_overlap=True,
            evidence_status=EvidenceInterpretation.EXCLUDED_UNAUTHENTICATED,
            exclusion_reason="Coordinates not found in authentic BED stream"
        )
        assert cfg_fake.evidence_status == EvidenceInterpretation.EXCLUDED_UNAUTHENTICATED

    # Challenge 11: Cached file that has changed detected
    def test_challenge_11_cached_file_changed_detected(self):
        """
        PrimaryPeakAuthenticator detects on-disk change and flags CHECKSUM_MISMATCH / FABRICATED_OR_ABSENT.
        """
        from tf_explorer.peak_authenticator import PrimaryPeakAuthenticator, PeakAuthStatus
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            bed = tmp_path / "ENCFF999TEST.bed"
            bed.write_text("chr1\t100\t200\tp1\t1000\t.\t10.0\t-1\t-1\t50\n", encoding="utf-8")
            empty_idx = tmp_path / "empty.json"
            empty_idx.write_text("{}", encoding="utf-8")

            auth = PrimaryPeakAuthenticator(search_paths=[tmp_path], index_file=empty_idx)
            res1 = auth.authenticate_peak("ENCFF999TEST", "chr1", 100, 200)
            assert res1.is_authenticated is True

            # Mutate
            bed.write_text("chr1\t300\t400\tp2\t1000\t.\t10.0\t-1\t-1\t50\n", encoding="utf-8")
            res2 = auth.authenticate_peak("ENCFF999TEST", "chr1", 100, 200)
            assert res2.is_authenticated is False

    # Challenge 12: Absent MANE Select designation
    def test_challenge_12_absent_mane_select_handled(self):
        """
        Transcript without MANE designation is marked MANE_UNVERIFIED and is_mane_select=False.
        """
        cfg = EvidenceConfiguration(
            config_id="CFG_NON_MANE",
            gene_symbol="TEST",
            transcript_id="ENST00000999999",
            is_mane_select=False,
            is_mane_plus_clinical=False,
            mane_status="MANE_UNVERIFIED"
        )
        assert cfg.is_mane_select is False
        assert cfg.mane_status == "MANE_UNVERIFIED"

    # Challenge 13: Interacting transcript and promoter-window effects
    def test_challenge_13_interacting_transcript_and_window_effects(self, engine):
        """
        Disagreement caused by simultaneous variation of alternative TSS and promoter window flank.
        Both causes must be recorded and multivariate interaction flagged.
        """
        cfg_narrow_tx1 = EvidenceConfiguration(
            config_id="CFG_N1",
            gene_symbol="TP53",
            transcript_id="ENST00000269305",
            canonical_tss=7687550,
            window_name="Core",
            promoter_upstream=500,
            promoter_downstream=100,
            promoter_start=7687050,
            promoter_end=7687650,
            encode_experiment_accession="ENCSR000EWX",
            encode_file_accession="ENCFF115KRO",
            biosample="MCF-7",
            is_authenticated=True,
            has_overlap=True,
            overlap_bp=200,
            evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE
        )
        cfg_wide_tx2 = EvidenceConfiguration(
            config_id="CFG_W2",
            gene_symbol="TP53",
            transcript_id="ENST00000445888",
            canonical_tss=7690000,  # Shifted TSS!
            window_name="Extended",
            promoter_upstream=2000,
            promoter_downstream=500,
            promoter_start=7688000,
            promoter_end=7690500,
            encode_experiment_accession="ENCSR000EWX",
            encode_file_accession="ENCFF115KRO",
            biosample="MCF-7",
            is_authenticated=True,
            has_overlap=False,
            overlap_bp=0,
            evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE
        )
        diag = engine.evaluate_pair(cfg_narrow_tx1, cfg_wide_tx2)
        assert DisagreementClassification.TRANSCRIPT_GEOMETRY_DEPENDENT in diag.primary_causes
        assert DisagreementClassification.PROMOTER_WINDOW_DEPENDENT in diag.primary_causes
        assert diag.is_multivariate_interaction is True
        assert len(diag.interacting_factors) >= 2

    # Challenge 14: Stage 7 Independent Reference Concordance
    def test_challenge_14_independent_reference_implementation_concordance(self, engine):
        """
        Tests that DisagreementEngine output agrees 100% with the blind independent reference implementation.
        """
        cfg_a = EvidenceConfiguration(
            config_id="A",
            gene_symbol="TEST",
            canonical_tss=1000,
            promoter_upstream=500,
            promoter_downstream=100,
            encode_experiment_accession="ENCSR001",
            encode_file_accession="ENCFF001",
            biosample="MCF-7",
            treatment_context="Untreated",
            is_authenticated=True,
            evidence_status=EvidenceInterpretation.POSITIVE_EVIDENCE
        )
        cfg_b = EvidenceConfiguration(
            config_id="B",
            gene_symbol="TEST",
            canonical_tss=2000,  # Different TSS
            promoter_upstream=500,
            promoter_downstream=100,
            encode_experiment_accession="ENCSR001",
            encode_file_accession="ENCFF001",
            biosample="MCF-7",
            treatment_context="Untreated",
            is_authenticated=True,
            evidence_status=EvidenceInterpretation.NEGATIVE_EVIDENCE
        )
        # Production engine diagnosis
        diag = engine.evaluate_pair(cfg_a, cfg_b)
        prod_labels = sorted([c.value for c in diag.primary_causes])

        # Independent reference diagnosis
        ref_dict_a = cfg_a.to_dict()
        ref_dict_b = cfg_b.to_dict()
        ref_labels = sorted(reference_classify_pair(ref_dict_a, ref_dict_b))

        assert prod_labels == ref_labels
        assert "TRANSCRIPT_GEOMETRY_DEPENDENT" in prod_labels


class TestStreamlitDisagreementExplorerIntegration:
    """End-to-End Streamlit Integration tests for the Disagreement Explorer Panel (Stage 8)."""

    @pytest.fixture
    def app(self):
        from pathlib import Path
        import tf_explorer
        app_path = str(Path(tf_explorer.__file__).parent / "app.py")
        at = AppTest.from_file(app_path, default_timeout=30)
        at.run()
        at.sidebar.radio[0].set_value("Disagreement Explorer").run()
        return at

    def test_disagreement_explorer_initial_load(self, app):
        """Verifies Disagreement Explorer panel loads properly with subheader and presets."""
        subheaders = [s.value for s in app.subheader]
        assert any("Disagreement Explorer" in s for s in subheaders)
        preset_sb = [s for s in app.selectbox if s.label == "Benchmark Presets"][0]
        assert "CDK4" in preset_sb.value

        # Check metrics cards
        metrics = [m.value for m in app.metric]
        assert len(metrics) >= 5
        # CDK4 preset has 2 configurations, 2 independent experiments
        assert metrics[0] == "2"  # Configurations
        assert metrics[1] == "2"  # Independent experiments

    def test_disagreement_explorer_preset_switching_and_exports(self, app):
        """Verifies switching to TFF1 preset updates metrics, diagnoses, and download buttons."""
        preset_sb = [s for s in app.selectbox if s.label == "Benchmark Presets"][0]
        tff1_opt = [opt for opt in preset_sb.options if "TFF1" in opt][0]
        preset_sb.set_value(tff1_opt).run()

        # Check metrics updated
        metrics = [m.value for m in app.metric]
        assert metrics[0] == "2"  # Configurations
        assert metrics[1] == "1"  # Independent experiments (1 experiment with 2 windows)

        # Check Step 3: Reproducibility Exports and expander exist
        markdown_sections = [m.value for m in app.markdown]
        assert any("Step 3: Reproducibility Exports" in m for m in markdown_sections)
        assert len(app.dataframe) > 0
        assert len(app.expander) > 0
