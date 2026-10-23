"""
Automated Regression and Unit Tests for Phase 3B Evidence Reliability Engine.
Tests EvidenceIntegrityEngine and RobustnessAnalyzer modules.
"""

import pytest
from tf_explorer.integrity import EvidenceIntegrityEngine, ValidationStatusCode
from tf_explorer.robustness import RobustnessAnalyzer, RobustnessClassification


class TestEvidenceIntegrityEngine:
    @pytest.fixture
    def engine(self):
        return EvidenceIntegrityEngine()

    def test_clean_validation(self, engine):
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX",
            expected_tf="E2F1",
            expected_biosample="MCF-7",
            expected_assembly="GRCh38",
            offline_only=False
        )
        assert res.is_valid is True
        assert res.status == ValidationStatusCode.VALIDATED
        assert "E2F1" in res.message

    def test_detects_tf_mismatch(self, engine):
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX",
            expected_tf="CTCF",
            expected_biosample="MCF-7",
            expected_assembly="GRCh38"
        )
        assert res.is_valid is False
        assert res.status == ValidationStatusCode.TF_TARGET_MISMATCH
        assert res.failure_code == "TF_TARGET_MISMATCH"

    def test_detects_file_experiment_mismatch(self, engine):
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000BSE",
            expected_tf="E2F1",
            expected_biosample="MCF-7",
            expected_assembly="GRCh38"
        )
        assert res.is_valid is False
        assert res.status == ValidationStatusCode.FILE_EXPERIMENT_MISMATCH

    def test_detects_biosample_mismatch(self, engine):
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX",
            expected_tf="E2F1",
            expected_biosample="K562",
            expected_assembly="GRCh38"
        )
        assert res.is_valid is False
        assert res.status == ValidationStatusCode.BIOSAMPLE_MISMATCH

    def test_detects_cross_experiment_biosample_mismatch(self, engine):
        """Verifies rejection when authentic MCF-7 CTCF dataset is queried against HCT116 expectation."""
        res = engine.verify(
            file_accession="ENCFF685HMV",
            experiment_accession="ENCSR000DMS",
            expected_tf="CTCF",
            expected_biosample="HCT116",
            expected_assembly="GRCh38"
        )
        assert res.is_valid is False
        assert res.status == ValidationStatusCode.BIOSAMPLE_MISMATCH


    def test_detects_assembly_mismatch(self, engine):
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX",
            expected_tf="E2F1",
            expected_biosample="MCF-7",
            expected_assembly="GRCh37"
        )
        assert res.is_valid is False
        assert res.status == ValidationStatusCode.ASSEMBLY_MISMATCH

    def test_detects_invalid_accession_format(self, engine):
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX_MCF7",
            expected_tf="E2F1",
            expected_biosample="MCF-7"
        )
        assert res.is_valid is False
        assert res.status == ValidationStatusCode.INVALID_ACCESSION_FORMAT

    def test_checksum_status_separation(self, engine, tmp_path):
        from tf_explorer.integrity import ChecksumStatusCode

        # 1. No file or expected MD5 -> FILE_NOT_DOWNLOADED
        res1 = engine.verify("ENCFF115KRO", "ENCSR000EWX", "E2F1", "MCF-7", "GRCh38")
        assert res1.checksum_status == ChecksumStatusCode.FILE_NOT_DOWNLOADED

        # 2. Remote checksum retrieved
        res2 = engine.verify(
            "ENCFF115KRO", "ENCSR000EWX", "E2F1", "MCF-7", "GRCh38",
            expected_md5sum="eeacedeabd905942553703034969b287"
        )
        assert res2.checksum_status == ChecksumStatusCode.REMOTE_CHECKSUM_RETRIEVED

        # 3. Local file hash mismatch
        dummy_file = tmp_path / "test.bed"
        dummy_file.write_text("chr1 100 200 peak1")
        res3 = engine.verify(
            "ENCFF115KRO", "ENCSR000EWX", "E2F1", "MCF-7", "GRCh38",
            local_file_path=str(dummy_file)
        )
        assert res3.is_valid is False
        assert res3.checksum_status == ChecksumStatusCode.CHECKSUM_MISMATCH

    def test_cache_schema_v2_attributes(self, engine):
        meta, src = engine.get_file_metadata("ENCFF115KRO")
        assert meta is not None
        # Must contain standard provenance attributes
        assert "file_accession" in meta
        assert "assembly" in meta
        assert meta["assembly"] == "GRCh38"


class TestRobustnessAnalyzer:
    @pytest.fixture
    def analyzer(self):
        return RobustnessAnalyzer()

    def test_strand_aware_coordinate_math(self, analyzer):
        # Positive strand: TSS=1000 [+], Up=500, Down=100 -> [500, 1100]
        p_pos = analyzer.calculate_promoter_coordinates(1000, "+", 500, 100)
        assert p_pos == (500, 1100)

        # Negative strand: TSS=58141400 [-], Up=1000, Down=200 -> [58141200, 58142400]
        p_neg = analyzer.calculate_promoter_coordinates(58141400, "-", 1000, 200)
        assert p_neg == (58141200, 58142400)

    def test_peak_overlap_math(self, analyzer):
        # Overlapping: peak [900, 1100] with promoter [1000, 2000] -> overlaps by 100 bp
        is_ov, bp = analyzer.check_peak_overlap(900, 1100, 1000, 2000)
        assert is_ov is True
        assert bp == 100

        # Non-overlapping abutment: peak [500, 1000] with promoter [1000, 2000] -> 0 bp
        is_ov, bp = analyzer.check_peak_overlap(500, 1000, 1000, 2000)
        assert is_ov is False
        assert bp == 0

    def test_transcript_sensitivity_classification(self, analyzer):
        txs = [
            {'id': 'ENST01', 'name': 'TX_PROX', 'tss': 1000, 'strand': '+', 'is_canonical': True},
            {'id': 'ENST02', 'name': 'TX_DIST', 'tss': 50000, 'strand': '+', 'is_canonical': False}
        ]
        exps = [{'experiment_id': 'EXP01', 'quality_tier': 'optimal', 'peaks': [{'start': 950, 'end': 1150}]}]
        rep = analyzer.evaluate("TEST_GENE", "TF1", "MCF-7", txs, exps, promoter_windows=[(500, 100)])
        assert rep.classification == RobustnessClassification.SENSITIVE_TO_TRANSCRIPT
        assert rep.transcript_support['ENST01'] is True
        assert rep.transcript_support['ENST02'] is False

    def test_experiment_sensitivity_classification(self, analyzer):
        txs = [{'id': 'ENST01', 'name': 'TX1', 'tss': 1000, 'strand': '+', 'is_canonical': True}]
        exps = [
            {'experiment_id': 'EXP01', 'quality_tier': 'optimal', 'peaks': [{'start': 950, 'end': 1150}]},
            {'experiment_id': 'EXP02', 'quality_tier': 'optimal', 'peaks': []}
        ]
        rep = analyzer.evaluate("TEST_GENE", "TF1", "MCF-7", txs, exps, promoter_windows=[(500, 100)])
        assert rep.classification == RobustnessClassification.SENSITIVE_TO_EXPERIMENT
        assert rep.experiment_support_fraction == "1/2"

    def test_window_parameter_sensitivity(self, analyzer):
        txs = [{'id': 'ENST01', 'name': 'TX1', 'tss': 2000, 'strand': '+', 'is_canonical': True}]
        # [-500, +100] -> [1500, 2100]; [-2000, +500] -> [0, 2500]
        # Peak at [200, 400] overlaps only the [-2000, +500] window
        exps = [{'experiment_id': 'EXP01', 'quality_tier': 'optimal', 'peaks': [{'start': 200, 'end': 400}]}]
        rep = analyzer.evaluate("TEST_GENE", "TF1", "MCF-7", txs, exps, promoter_windows=[(500, 100), (2000, 500)])
        assert rep.classification == RobustnessClassification.SENSITIVE_TO_PARAMETERS

    def test_insufficient_data_classification(self, analyzer):
        txs = [{'id': 'ENST01', 'name': 'TX1', 'tss': 1000, 'strand': '+', 'is_canonical': True}]
        rep = analyzer.evaluate("TEST_GENE", "TF1", "MCF-7", txs, experiment_peak_sets=[])
        assert rep.classification == RobustnessClassification.INSUFFICIENT_DATA

    def test_real_biology_cdk4_grch38_transcript_sensitivity(self, analyzer):
        """Regression test: Authentic GRCh38 CDK4 transcript sensitivity (Negative strand)."""
        txs = [
            {"id": "ENST00000257904.11", "name": "CDK4-201 (Canonical)", "tss": 57752310, "strand": "-", "is_canonical": True},
            {"id": "ENST00000552713.5", "name": "CDK4-204 (Alternative)", "tss": 57749795, "strand": "-", "is_canonical": False}
        ]
        exps = [{
            "experiment_id": "ENCSR000EWX",
            "quality_tier": "optimal",
            "peaks": [{"start": 57752169, "end": 57752512, "id": "cdk4_pk1"}]
        }]
        rep = analyzer.evaluate("CDK4", "E2F1", "MCF-7", txs, exps, promoter_windows=[(500, 100), (1000, 200)])
        assert rep.classification == RobustnessClassification.SENSITIVE_TO_TRANSCRIPT
        assert rep.transcript_support["ENST00000257904.11"] is True
        assert rep.transcript_support["ENST00000552713.5"] is False

    def test_real_biology_tp53_grch38_parameter_sensitivity(self, analyzer):
        """Regression test: Authentic GRCh38 TP53 parameter sensitivity (CTCF upstream peak at -704 bp)."""
        txs = [
            {"id": "ENST00000269305.9", "name": "TP53-201 (Canonical)", "tss": 7687490, "strand": "-", "is_canonical": True}
        ]
        exps = [{
            "experiment_id": "ENCSR000DMS",
            "quality_tier": "optimal",
            "peaks": [{"start": 7688194, "end": 7688248, "id": "tp53_ctcf_pk1"}]
        }]
        # Window [-500, +100] -> [7687390, 7687990] (peak 7688194 is outside)
        # Window [-1000, +200] -> [7687290, 7688490] (peak 7688194 is inside)
        rep = analyzer.evaluate("TP53", "CTCF", "MCF-7", txs, exps, promoter_windows=[(500, 100), (1000, 200)])
        assert rep.classification == RobustnessClassification.SENSITIVE_TO_PARAMETERS
        assert rep.window_support["-500/+100"] is False
        assert rep.window_support["-1000/+200"] is True

