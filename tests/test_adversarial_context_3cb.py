"""
Phase 3C-B: Stage 2 — Experimental-Context Adversarial Tests (Scenarios A through G).
Tests strict unmodified, strict untreated, explicit acceptance, unspecified requirements,
incomplete metadata, exploratory override, and replicate discordance.
"""

import pytest
from tf_explorer.integrity import (
    EvidenceIntegrityEngine,
    ValidationStatusCode,
    ExperimentalContextStatus,
    ValidationResult
)


class TestExperimentalContextAdversarialScenarios:
    @pytest.fixture
    def engine(self):
        return EvidenceIntegrityEngine()

    def test_scenario_a_strict_unmodified_context(self, engine):
        """Scenario A: Authentic ENCSR000EWX has N-terminal HA tag; strict unmodified must trigger CONTEXT_MISMATCH."""
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX",
            expected_tf="E2F1",
            expected_biosample="MCF-7",
            expected_unmodified=True,
            allow_exploratory_context=False
        )
        assert res.is_valid is True  # Authenticity valid
        assert res.status == ValidationStatusCode.VALIDATED
        assert res.context_status == ExperimentalContextStatus.CONTEXT_MISMATCH
        assert res.context_compatible is False
        assert any("genetically modified" in w.lower() or "ha" in w.lower() for w in res.context_warnings)

    def test_scenario_b_strict_untreated_context(self, engine):
        """Scenario B: Authentic ENCSR000DMS has 100 nM 17β-estradiol; strict untreated must trigger CONTEXT_MISMATCH."""
        res = engine.verify(
            file_accession="ENCFF685HMV",
            experiment_accession="ENCSR000DMS",
            expected_tf="CTCF",
            expected_biosample="MCF-7",
            expected_untreated=True,
            allow_exploratory_context=False
        )
        assert res.is_valid is True  # Authenticity valid
        assert res.status == ValidationStatusCode.VALIDATED
        assert res.context_status == ExperimentalContextStatus.CONTEXT_MISMATCH
        assert res.context_compatible is False
        assert any("treated" in w.lower() or "estradiol" in w.lower() for w in res.context_warnings)

    def test_scenario_c_genetically_modified_context_explicitly_accepted(self, engine):
        """Scenario C: Requested tagged ESR1 matches CRISPR eGFP tag in ENCSR463GOT -> CONTEXT_MATCH."""
        res = engine.verify(
            file_accession="ENCFF185UUT",
            experiment_accession="ENCSR463GOT",
            expected_tf="ESR1",
            expected_biosample="MCF-7",
            genetic_modification_requirement="eGFP",
            allow_exploratory_context=False
        )
        assert res.is_valid is True
        assert res.context_status == ExperimentalContextStatus.CONTEXT_MATCH
        assert res.context_compatible is True

    def test_scenario_d_context_requirement_unspecified(self, engine):
        """Scenario D: When requirements are unspecified, engine yields CONTEXT_UNSPECIFIED, not false equivalence."""
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX",
            expected_tf="E2F1",
            expected_biosample="MCF-7",
            treatment_requirement=None,
            genetic_modification_requirement=None
        )
        assert res.is_valid is True
        assert res.context_status == ExperimentalContextStatus.CONTEXT_UNSPECIFIED
        # Engine attaches informative notice but does not claim biological equivalence was tested
        assert any("Note:" in w for w in res.context_warnings)

    def test_scenario_e_incomplete_metadata(self):
        """Scenario E: When experiment record lacks treatment and modification fields, return CONTEXT_METADATA_INCOMPLETE."""
        incomplete_snapshot = {
            "encode_files": {
                "ENCFF000INC": {
                    "experiment": "ENCSR000INC",
                    "target": "MYC",
                    "biosample": "K562",
                    "assembly": "GRCh38",
                    "format": "bed",
                    "output_type": "narrowPeak",
                    "status": "released"
                }
            },
            "encode_experiments": {
                "ENCSR000INC": {
                    "target": "MYC",
                    "biosample": "K562",
                    "status": "released"
                    # Neither is_treated, treatments, is_genetically_modified, nor genetic_modifications
                }
            }
        }
        engine = EvidenceIntegrityEngine(memory_snapshot=incomplete_snapshot)
        res = engine.verify(
            file_accession="ENCFF000INC",
            experiment_accession="ENCSR000INC",
            expected_tf="MYC",
            expected_biosample="K562",
            expected_untreated=True,
            expected_unmodified=True,
            allow_exploratory_context=False
        )
        assert res.is_valid is True
        assert res.context_status == ExperimentalContextStatus.CONTEXT_METADATA_INCOMPLETE
        assert res.context_compatible is False
        assert any("lacks authoritative" in w.lower() for w in res.context_warnings)

    def test_scenario_f_exploratory_override(self, engine):
        """Scenario F: Context mismatch with allow_exploratory_context=True yields CONTEXT_ACCEPTED_WITH_LIMITATIONS."""
        res = engine.verify(
            file_accession="ENCFF115KRO",
            experiment_accession="ENCSR000EWX",
            expected_tf="E2F1",
            expected_biosample="MCF-7",
            expected_unmodified=True,
            allow_exploratory_context=True
        )
        assert res.is_valid is True
        assert res.context_status == ExperimentalContextStatus.CONTEXT_ACCEPTED_WITH_LIMITATIONS
        assert res.context_compatible is True
        assert res.context_accepted is True
        # Preserves original mismatch warning
        assert any("genetically modified" in w.lower() for w in res.context_warnings)
        # Preserves exploratory qualification warning
        assert any("EXPLORATORY CONTEXT QUALIFICATION" in w for w in res.context_warnings)

    def test_scenario_g_replicate_discordance(self):
        """Scenario G: Discordant replicates (Rep 1 untreated, Rep 2 treated) must be explicitly flagged."""
        discordant_snapshot = {
            "encode_files": {
                "ENCFF000DIS": {
                    "experiment": "ENCSR000DIS",
                    "target": "JUN",
                    "biosample": "HepG2",
                    "assembly": "GRCh38",
                    "format": "bed",
                    "output_type": "narrowPeak",
                    "status": "released"
                }
            },
            "encode_experiments": {
                "ENCSR000DIS": {
                    "target": "JUN",
                    "biosample": "HepG2",
                    "status": "released",
                    "is_treated": True,
                    "treatments": ["ethanol"],
                    "is_genetically_modified": False,
                    "genetic_modifications": [],
                    "replicates": [
                        {
                            "biological_replicate_number": 1,
                            "treatments": [],
                            "genetic_modifications": [],
                            "is_treated": False,
                            "is_genetically_modified": False
                        },
                        {
                            "biological_replicate_number": 2,
                            "treatments": ["ethanol"],
                            "genetic_modifications": [],
                            "is_treated": True,
                            "is_genetically_modified": False
                        }
                    ]
                }
            }
        }
        engine = EvidenceIntegrityEngine(memory_snapshot=discordant_snapshot)
        res = engine.verify(
            file_accession="ENCFF000DIS",
            experiment_accession="ENCSR000DIS",
            expected_tf="JUN",
            expected_biosample="HepG2"
        )
        assert res.is_valid is True
        assert res.biological_context["has_replicate_discordance"] is True
        assert len(res.biological_context["replicates"]) == 2
        assert any("REPLICATE CONTEXT DISCORDANCE" in w for w in res.context_warnings)
