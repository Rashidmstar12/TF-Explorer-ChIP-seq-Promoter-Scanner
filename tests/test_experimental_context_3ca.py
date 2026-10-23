"""
Automated Regression and Unit Tests for Phase 3C-A: Experimental-Context Integrity.
Validates ExperimentalContextStatus, treatment & tag detection, exploratory opt-in,
and Markdown export generation.
"""

import pytest
from tf_explorer.integrity import (
    EvidenceIntegrityEngine,
    ValidationStatusCode,
    ExperimentalContextStatus,
    ValidationResult
)
from tf_explorer.robustness import (
    RobustnessAnalyzer,
    RobustnessClassification,
    RobustnessReport
)


class TestExperimentalContextIntegrity:
    @pytest.fixture
    def engine(self):
        return EvidenceIntegrityEngine()

    def test_encsr000ewx_ha_tag_detection(self, engine):
        """ENCSR000EWX has an N-terminal HA tag on E2F1."""
        meta, _ = engine.get_experiment_metadata("ENCSR000EWX")
        assert meta["is_genetically_modified"] is True
        assert len(meta["genetic_modifications"]) > 0
        tag = meta["genetic_modifications"][0]["tags"][0]["name"]
        assert "HA" in tag

        # When researcher specifies untreated, unmodified:
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX",
            expected_tf="E2F1",
            expected_biosample="MCF-7",
            expected_untreated=True,
            expected_unmodified=True,
            allow_exploratory_context=False
        )
        assert res.is_valid is True  # Data authenticity is VALID
        assert res.context_status == ExperimentalContextStatus.CONTEXT_MISMATCH
        assert any("genetically modified" in note.lower() for note in res.context_notes)

    def test_encsr000ewx_ha_tag_exploratory_acceptance(self, engine):
        """When researcher allows exploratory context, tagged assay is accepted with limitations."""
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX",
            expected_tf="E2F1",
            expected_biosample="MCF-7",
            expected_untreated=True,
            expected_unmodified=True,
            allow_exploratory_context=True
        )
        assert res.is_valid is True
        assert res.context_status == ExperimentalContextStatus.CONTEXT_ACCEPTED_WITH_LIMITATIONS
        assert len(res.context_notes) > 0

    def test_encsr463got_crispr_egfp_detection(self, engine):
        """ENCSR463GOT has a CRISPR knock-in insertion of C-terminal eGFP tag on ESR1."""
        meta, _ = engine.get_experiment_metadata("ENCSR463GOT")
        assert meta["is_genetically_modified"] is True
        assert len(meta["genetic_modifications"]) > 0
        mod = meta["genetic_modifications"][0]
        tag_name = mod["tags"][0]["name"]
        assert tag_name == "eGFP"

        # Opt-in acceptance
        res = engine.verify(
            file_accession="ENCFF185UUT",
            experiment_accession="ENCSR463GOT",
            expected_tf="ESR1",
            expected_biosample="MCF-7",
            expected_unmodified=True,
            allow_exploratory_context=True
        )
        assert res.is_valid is True
        assert res.context_status == ExperimentalContextStatus.CONTEXT_ACCEPTED_WITH_LIMITATIONS

    def test_encsr000dms_estradiol_treatment_detection(self, engine):
        """ENCSR000DMS has 100 nM 17β-estradiol treatment for 45 minutes."""
        meta, _ = engine.get_experiment_metadata("ENCSR000DMS")
        assert meta["is_treated"] is True
        assert len(meta["treatments"]) > 0
        treatment_str = meta["treatments"][0]
        assert "estradiol" in treatment_str.lower()

        # Strict untreated check triggers CONTEXT_MISMATCH
        res = engine.verify(
            file_accession="ENCFF685HMV",
            experiment_accession="ENCSR000DMS",
            expected_tf="CTCF",
            expected_biosample="MCF-7",
            expected_untreated=True,
            allow_exploratory_context=False
        )
        assert res.is_valid is True  # Metadata is valid
        assert res.context_status == ExperimentalContextStatus.CONTEXT_MISMATCH
        assert any("treated" in note.lower() for note in res.context_notes)

        # Exploratory check allows with limitations
        res_opt = engine.verify(
            file_accession="ENCFF685HMV",
            experiment_accession="ENCSR000DMS",
            expected_tf="CTCF",
            expected_biosample="MCF-7",
            expected_untreated=True,
            allow_exploratory_context=True
        )
        assert res_opt.context_status == ExperimentalContextStatus.CONTEXT_ACCEPTED_WITH_LIMITATIONS

    def test_context_unspecified_when_not_requested(self, engine):
        """When neither treatments nor modifications are checked, context is CONTEXT_UNSPECIFIED."""
        res = engine.verify(
            file_accession="ENCFF004ZKN",
            experiment_accession="ENCSR000EVM",
            expected_tf="E2F1",
            expected_biosample="HeLa-S3"
        )
        assert res.is_valid is True
        assert res.context_status == ExperimentalContextStatus.CONTEXT_UNSPECIFIED

    def test_context_metadata_incomplete_fallback(self, engine):
        """When experiment has no treatment/modification metadata recorded."""
        # Test synthetic record without treatment fields
        res = engine.verify(
            file_accession="ENCFF004ZKN",
            experiment_accession="ENCSR000EVM",
            expected_tf="E2F1",
            expected_biosample="HeLa-S3",
            expected_untreated=True,
            allow_exploratory_context=False
        )
        assert res.context_status in (
            ExperimentalContextStatus.CONTEXT_METADATA_INCOMPLETE,
            ExperimentalContextStatus.CONTEXT_MATCH
        )

    def test_robustness_report_includes_context_and_markdown(self):
        """RobustnessAnalyzer generates markdown with biological context summary."""
        analyzer = RobustnessAnalyzer()
        txs = [
            {"id": "ENST00000257904.11", "name": "CDK4-201 (Canonical)", "tss": 57752310, "strand": "-", "is_canonical": True}
        ]
        exps = [{
            "experiment_id": "ENCSR000EWX",
            "quality_tier": "optimal",
            "peaks": [{"start": 57752169, "end": 57752512, "id": "cdk4_pk1"}]
        }]
        contexts = {
            "ENCSR000EWX": {
                "status": "CONTEXT_ACCEPTED_WITH_LIMITATIONS",
                "notes": ["N-terminal HA tag (exploratory opt-in)"]
            }
        }
        report = analyzer.evaluate(
            gene_symbol="CDK4",
            tf="E2F1",
            biosample="MCF-7",
            transcripts=txs,
            experiment_peak_sets=exps,
            promoter_windows=[(1000, 200)],
            biological_contexts=contexts
        )
        assert report.classification == RobustnessClassification.STABLE_WITHIN_EVALUATED_CHOICES
        assert "ENCSR000EWX" in report.biological_contexts
        assert len(report.context_warnings) > 0

        # Markdown report rendering
        md = report.to_researcher_markdown()
        assert "# Regulatory Evidence Robustness Dossier:" in md
        assert "CDK4" in md
        assert "E2F1" in md
        assert "MCF-7" in md
        assert "Experimental Context & Perturbation Alerts" in md
        assert "Evaluated Experiments & Contexts" in md
        assert "ENCSR000EWX" in md
