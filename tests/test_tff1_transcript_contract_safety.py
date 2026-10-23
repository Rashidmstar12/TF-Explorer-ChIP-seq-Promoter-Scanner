"""
Targeted Positive and Negative Regression Tests for TFF1 Transcript Contract Safety.
Ensures strict biological authentication for canonical MANE Select ENST00000291527.3,
fail-closed rejection of legacy .9 and unknown versions, and non-circular derivation
across both the high-level Factory API and direct EvidenceContractValidator layers.

Distinguishes explicitly between:
- Scenario 1: Offline synthetic geometry tests (has_synthetic_configuration=True)
- Scenario 2: Missing source data handling (SOURCE_UNAVAILABLE)
- Scenario 3: Authentic real-data integration against primary ENCODE BED (requires physical file)
"""

import os
from pathlib import Path
from typing import Optional
import pytest

from tf_explorer.disagreement_factory import (
    build_validated_configuration,
    find_source_bed_file,
    create_synthetic_configuration_for_testing
)
from tf_explorer.disagreement_provenance import resolve_transcript_version
from tf_explorer.disagreement_taxonomy import (
    EvidenceInterpretation,
    DisagreementClassification,
    DimensionType,
    EvidenceConfiguration
)
from tf_explorer.disagreement_engine import DisagreementEngine
from tf_explorer.evidence_contract import (
    EvidenceContractValidator,
    EvidenceRecord,
    ValidationStatus
)


def _get_tff1_source_bed() -> Optional[Path]:
    """Locate the authentic historical ENCFF185UUT peak dataset if available on this system."""
    env_path = os.environ.get("TF_EXPLORER_TFF1_BED_PATH")
    if env_path and Path(env_path).is_file():
        return Path(env_path)

    candidates = [
        Path.cwd() / "reference_data" / "peaks" / "ENCFF185UUT.bed.gz",
        Path(__file__).resolve().parent.parent / "reference_data" / "peaks" / "ENCFF185UUT.bed.gz",
    ]
    for c in candidates:
        if c.is_file():
            return c

    bf = find_source_bed_file("ENCFF185UUT")
    if bf and bf.is_file():
        return bf

    return None


class TestTFF1TranscriptContractSafety:
    """Verifies that TFF1 transcript version checking is scientifically defensible and fail-closed."""

    def test_01_canonical_mane_select_v3_authenticates_cleanly(self):
        """Canonical MANE Select transcript ENST00000291527.3 must authenticate and detect primary peak."""
        resolved, note = resolve_transcript_version("TFF1", "ENST00000291527.3")
        assert resolved == "ENST00000291527.3"
        assert note is None

        bed_file = _get_tff1_source_bed()
        cfg = build_validated_configuration(
            gene_symbol="TFF1",
            transcript_id="ENST00000291527.3",
            canonical_tss=42366535,
            strand="-",
            tf_target="ESR1",
            encode_experiment_accession="ENCSR463GOT",
            encode_file_accession="ENCFF185UUT",
            biosample="MCF-7",
            window_name="Standard_Proximal",
            window_start_offset=-2000,
            window_end_offset=500,
            chromosome="chr21",
            custom_bed_path=bed_file
        )
        assert cfg.transcript_id == "ENST00000291527.3"
        if bed_file and bed_file.exists():
            assert cfg.is_authenticated is True
            assert cfg.auth_status == "AUTHENTICATED_PRIMARY_BED"
            assert cfg.evidence_status == EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION
            assert cfg.has_overlap is True
            assert cfg.overlap_bp == 691
            assert cfg.exclusion_reason is None
        else:
            # In isolated CI runner without external data, verify clean fail-closed source unavailability
            assert cfg.is_authenticated is False
            assert cfg.auth_status == "SOURCE_UNAVAILABLE"
            assert cfg.evidence_status == EvidenceInterpretation.SOURCE_UNAVAILABLE

    def test_02_unversioned_resolves_to_canonical_v3_with_audit_note(self):
        """Unversioned identifier ENST00000291527 must resolve directly to canonical .3 with audit note."""
        resolved, note = resolve_transcript_version("TFF1", "ENST00000291527")
        assert resolved == "ENST00000291527.3"
        assert note == "Resolved unversioned transcript 'ENST00000291527' to authoritative version 'ENST00000291527.3'"

        bed_file = _get_tff1_source_bed()
        cfg = build_validated_configuration(
            gene_symbol="TFF1",
            transcript_id="ENST00000291527",
            canonical_tss=42366535,
            strand="-",
            tf_target="ESR1",
            encode_experiment_accession="ENCSR463GOT",
            encode_file_accession="ENCFF185UUT",
            biosample="MCF-7",
            window_name="Standard_Proximal",
            window_start_offset=-2000,
            window_end_offset=500,
            chromosome="chr21",
            custom_bed_path=bed_file
        )
        assert cfg.transcript_id == "ENST00000291527.3"
        assert any("Resolved unversioned transcript" in n for n in cfg.notes)
        if bed_file and bed_file.exists():
            assert cfg.is_authenticated is True
            assert cfg.auth_status == "AUTHENTICATED_PRIMARY_BED"
            assert cfg.evidence_status == EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION
        else:
            assert cfg.is_authenticated is False
            assert cfg.auth_status == "SOURCE_UNAVAILABLE"
            assert cfg.evidence_status == EvidenceInterpretation.SOURCE_UNAVAILABLE

    def test_03_legacy_v9_fails_closed_without_biological_authentication(self):
        """Legacy fixture version .9 must NOT obtain biological authentication merely for backward compatibility."""
        resolved, note = resolve_transcript_version("TFF1", "ENST00000291527.9")
        assert resolved is None
        assert "requested '9', authoritative is '3'" in note

        bed_file = _get_tff1_source_bed()
        cfg = build_validated_configuration(
            gene_symbol="TFF1",
            transcript_id="ENST00000291527.9",
            canonical_tss=42366535,
            strand="-",
            tf_target="ESR1",
            encode_experiment_accession="ENCSR463GOT",
            encode_file_accession="ENCFF185UUT",
            biosample="MCF-7",
            window_name="Standard_Proximal",
            window_start_offset=-2000,
            window_end_offset=500,
            chromosome="chr21",
            custom_bed_path=bed_file
        )
        # Preflight fails closed regardless of whether BED file exists
        assert cfg.is_authenticated is False
        assert cfg.auth_status == "GENE_TRANSCRIPT_MISMATCH"
        assert cfg.evidence_status == EvidenceInterpretation.SOURCE_NOT_VERIFIED
        assert "Transcript version mismatch" in cfg.exclusion_reason
        assert "requested '9', authoritative is '3'" in cfg.exclusion_reason
        assert cfg.transcript_id == "ENST00000291527.9"

    def test_04_arbitrary_unknown_versions_fail_closed(self):
        """Arbitrary version suffixes must fail closed fail-safe."""
        bed_file = _get_tff1_source_bed()
        for bad_ver in ["ENST00000291527.999", "ENST00000291527.1", "ENST00000291527.bad"]:
            resolved, note = resolve_transcript_version("TFF1", bad_ver)
            assert resolved is None
            assert "authoritative is '3'" in note

            cfg = build_validated_configuration(
                gene_symbol="TFF1",
                transcript_id=bad_ver,
                canonical_tss=42366535,
                strand="-",
                tf_target="ESR1",
                encode_experiment_accession="ENCSR463GOT",
                encode_file_accession="ENCFF185UUT",
                biosample="MCF-7",
                window_name="Standard",
                window_start_offset=-2000,
                window_end_offset=500,
                chromosome="chr21",
                custom_bed_path=bed_file
            )
            assert cfg.is_authenticated is False
            assert cfg.auth_status == "GENE_TRANSCRIPT_MISMATCH"
            assert cfg.evidence_status == EvidenceInterpretation.SOURCE_NOT_VERIFIED
            assert "Transcript version mismatch" in cfg.exclusion_reason
            assert cfg.transcript_id == bad_ver

    def test_05_offline_synthetic_geometry_scenario(self):
        """Scenario 1: Offline synthetic geometry test verifies classifier without empirical biological authentication."""
        cfg_std = create_synthetic_configuration_for_testing(
            config_id="TFF1_Standard_Window_Synthetic",
            gene_symbol="TFF1",
            transcript_id="ENST00000291527.3",
            canonical_tss=42366535,
            strand="-",
            window_name="Standard_Proximal",
            window_start_offset=-2000,
            window_end_offset=500,
            tf_target="ESR1",
            biosample="MCF-7",
            encode_experiment_accession="ENCSR463GOT",
            encode_file_accession="ENCFF185UUT",
            has_overlap=True,
            overlap_bp=691,
            evidence_status=EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION
        )
        cfg_dist = create_synthetic_configuration_for_testing(
            config_id="TFF1_Distal_Upstream_Synthetic",
            gene_symbol="TFF1",
            transcript_id="ENST00000291527.3",
            canonical_tss=42366535,
            strand="-",
            window_name="Distal_Upstream",
            window_start_offset=-3000,
            window_end_offset=-1000,
            tf_target="ESR1",
            biosample="MCF-7",
            encode_experiment_accession="ENCSR463GOT",
            encode_file_accession="ENCFF185UUT",
            has_overlap=False,
            overlap_bp=0,
            evidence_status=EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED
        )

        assert cfg_std.is_synthetic is True
        assert cfg_std.is_verified_biological_observation() is False
        assert cfg_dist.is_synthetic is True
        assert cfg_dist.is_verified_biological_observation() is False

        engine = DisagreementEngine()
        diag = engine.evaluate_pair(cfg_std, cfg_dist)

        assert diag.has_synthetic_configuration is True
        assert diag.primary_causes == [DisagreementClassification.PROMOTER_WINDOW_DEPENDENT]
        assert diag.differing_dimensions == [DimensionType.PROMOTER_WINDOW]
        assert "Promoter window boundary parameter changes" in diag.scientific_explanation

    def test_06_missing_real_source_scenario(self):
        """Scenario 2: Calling factory without physical BED produces fail-closed SOURCE_UNAVAILABLE."""
        missing_bed = Path("/nonexistent/data/ENCFF185UUT.bed.gz")
        cfg_std = build_validated_configuration(
            gene_symbol="TFF1",
            transcript_id="ENST00000291527.3",
            canonical_tss=42366535,
            strand="-",
            tf_target="ESR1",
            encode_experiment_accession="ENCSR463GOT",
            encode_file_accession="ENCFF185UUT",
            biosample="MCF-7",
            window_name="Standard_Proximal",
            window_start_offset=-2000,
            window_end_offset=500,
            chromosome="chr21",
            custom_bed_path=missing_bed
        )
        cfg_dist = build_validated_configuration(
            gene_symbol="TFF1",
            transcript_id="ENST00000291527.3",
            canonical_tss=42366535,
            strand="-",
            tf_target="ESR1",
            encode_experiment_accession="ENCSR463GOT",
            encode_file_accession="ENCFF185UUT",
            biosample="MCF-7",
            window_name="Distal_Upstream",
            window_start_offset=-3000,
            window_end_offset=-1000,
            chromosome="chr21",
            custom_bed_path=missing_bed
        )

        assert cfg_std.is_authenticated is False
        assert cfg_std.auth_status == "SOURCE_UNAVAILABLE"
        assert cfg_std.evidence_status == EvidenceInterpretation.SOURCE_UNAVAILABLE
        assert cfg_std.is_verified_biological_observation() is False

        assert cfg_dist.is_authenticated is False
        assert cfg_dist.auth_status == "SOURCE_UNAVAILABLE"
        assert cfg_dist.evidence_status == EvidenceInterpretation.SOURCE_UNAVAILABLE
        assert cfg_dist.is_verified_biological_observation() is False

        engine = DisagreementEngine()
        diag = engine.evaluate_pair(cfg_std, cfg_dist)
        assert diag.evidence_a == EvidenceInterpretation.SOURCE_UNAVAILABLE
        assert diag.evidence_b == EvidenceInterpretation.SOURCE_UNAVAILABLE
        assert "Comparison cannot establish biological or geometric concordance" in diag.scientific_explanation

    def test_07_authentic_real_data_integration_scenario(self):
        """Scenario 3: Authentic real-data integration executing only when verified primary BED is available."""
        bed_file = _get_tff1_source_bed()
        if bed_file is None or not bed_file.is_file():
            pytest.skip(
                "Authentic primary BED dataset ENCFF185UUT not found on local disk. "
                "Real-data biological integration test skipped in offline / clean CI environment."
            )

        cfg_std = build_validated_configuration(
            gene_symbol="TFF1",
            transcript_id="ENST00000291527.3",
            canonical_tss=42366535,
            strand="-",
            tf_target="ESR1",
            encode_experiment_accession="ENCSR463GOT",
            encode_file_accession="ENCFF185UUT",
            biosample="MCF-7",
            window_name="Standard_Proximal",
            window_start_offset=-2000,
            window_end_offset=500,
            chromosome="chr21",
            custom_bed_path=bed_file
        )
        cfg_dist = build_validated_configuration(
            gene_symbol="TFF1",
            transcript_id="ENST00000291527.3",
            canonical_tss=42366535,
            strand="-",
            tf_target="ESR1",
            encode_experiment_accession="ENCSR463GOT",
            encode_file_accession="ENCFF185UUT",
            biosample="MCF-7",
            window_name="Distal_Upstream",
            window_start_offset=-3000,
            window_end_offset=-1000,
            chromosome="chr21",
            custom_bed_path=bed_file
        )

        assert cfg_std.is_authenticated is True
        assert cfg_std.auth_status == "AUTHENTICATED_PRIMARY_BED"
        assert cfg_std.evidence_status == EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION
        assert cfg_std.has_overlap is True
        assert cfg_std.overlap_bp == 691
        assert cfg_std.is_verified_biological_observation() is True

        assert cfg_dist.is_authenticated is True
        assert cfg_dist.auth_status == "AUTHENTICATED_PRIMARY_BED"
        assert cfg_dist.evidence_status == EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED
        assert cfg_dist.has_overlap is False
        assert cfg_dist.overlap_bp == 0
        assert cfg_dist.is_verified_biological_observation() is True

        engine = DisagreementEngine()
        diag = engine.evaluate_pair(cfg_std, cfg_dist)

        assert diag.has_synthetic_configuration is False
        assert diag.primary_causes == [DisagreementClassification.PROMOTER_WINDOW_DEPENDENT]
        assert diag.differing_dimensions == [DimensionType.PROMOTER_WINDOW]
        assert "Promoter window boundary parameter changes" in diag.scientific_explanation

    def test_08_direct_evidence_contract_validator_rejects_tff1_v9(self):
        """Direct EvidenceContractValidator must fail closed on TFF1 version 9 without bypass."""
        validator = EvidenceContractValidator()
        rec_v9 = EvidenceRecord(
            gene_symbol="TFF1",
            ncbi_gene_id="7031",
            ensembl_gene_id="ENSG00000160182",
            transcript_id="ENST00000291527",
            transcript_version="9",
            genome_assembly="GRCh38",
            chromosome="chr21",
            strand="-",
            canonical_tss=42366535,
            promoter_upstream=2000,
            promoter_downstream=500,
            promoter_start=42366035,
            promoter_end=42368535,
            encode_file_accession="ENCFF185UUT",
            encode_experiment_accession="ENCSR463GOT",
            experiment_target="ESR1",
            biosample="MCF-7",
            treatment_context="Untreated",
            genetic_modification_context="eGFP-tagged (C-terminal)",
            source_peak_start=42366395,
            source_peak_end=42367086,
            source_peak_signal=2601.5,
            source_peak_summit=42366729,
            overlap_bp=691,
            overlap_percentage=100.0,
            evidence_source_url="https://www.encodeproject.org/files/ENCFF185UUT/",
            retrieval_date="2026-10-23",
            file_md5="464964b5a099bdd901dfab66de064b21",
            peak_chromosome="chr21"
        )
        res = validator.validate(rec_v9)
        assert res.is_valid is False
        assert res.status == ValidationStatus.FAILED
        assert any(
            "Gene 'TFF1' transcript version mismatch: Authoritative MANE Select version is '3', but record claims '9'."
            in failure for failure in res.failures
        )

    def test_09_direct_evidence_contract_validator_accepts_tff1_v3(self):
        """Direct EvidenceContractValidator must accept canonical TFF1 version 3."""
        bed_file = _get_tff1_source_bed()
        if bed_file and bed_file.exists():
            from tf_explorer.peak_authenticator import PrimaryPeakAuthenticator
            auth = PrimaryPeakAuthenticator(search_paths=[bed_file.parent])
            validator = EvidenceContractValidator(peak_authenticator=auth)
        else:
            validator = EvidenceContractValidator()

        rec_v3 = EvidenceRecord(
            gene_symbol="TFF1",
            ncbi_gene_id="7031",
            ensembl_gene_id="ENSG00000160182",
            transcript_id="ENST00000291527",
            transcript_version="3",
            genome_assembly="GRCh38",
            chromosome="chr21",
            strand="-",
            canonical_tss=42366535,
            promoter_upstream=2000,
            promoter_downstream=500,
            promoter_start=42366035,
            promoter_end=42368535,
            encode_file_accession="ENCFF185UUT",
            encode_experiment_accession="ENCSR463GOT",
            experiment_target="ESR1",
            biosample="MCF-7",
            treatment_context="Untreated",
            genetic_modification_context="eGFP-tagged (C-terminal)",
            source_peak_start=42366395,
            source_peak_end=42367086,
            source_peak_signal=2601.5,
            source_peak_summit=42366729,
            overlap_bp=691,
            overlap_percentage=100.0,
            evidence_source_url="https://www.encodeproject.org/files/ENCFF185UUT/",
            retrieval_date="2026-10-23",
            file_md5="464964b5a099bdd901dfab66de064b21",
            peak_chromosome="chr21"
        )
        res = validator.validate(rec_v3)
        assert len(res.failures) == 0
        assert not any("transcript version mismatch" in f for f in res.failures)
        if bed_file and bed_file.exists():
            assert res.is_valid is True
            assert res.status == ValidationStatus.VERIFIED_PASS

    def test_10_direct_evidence_contract_validator_other_genes_retained(self):
        """Direct EvidenceContractValidator must retain correct version validation on other genes."""
        validator = EvidenceContractValidator()
        # CDK4 canonical is version 11
        rec_cdk4_valid = EvidenceRecord(
            gene_symbol="CDK4",
            ncbi_gene_id="1019",
            ensembl_gene_id="ENSG00000135446",
            transcript_id="ENST00000257904",
            transcript_version="11",
            genome_assembly="GRCh38",
            chromosome="chr12",
            strand="-",
            canonical_tss=57752310,
            promoter_upstream=2000,
            promoter_downstream=500,
            promoter_start=57751810,
            promoter_end=57754310,
            encode_file_accession="ENCFF000TEST",
            encode_experiment_accession="ENCSR000TEST",
            experiment_target="CDK4",
            biosample="K562",
            treatment_context="Untreated",
            genetic_modification_context="None",
            source_peak_start=57752000,
            source_peak_end=57752500,
            source_peak_signal=50.0,
            source_peak_summit=57752250,
            overlap_bp=500,
            overlap_percentage=100.0,
            evidence_source_url="https://www.encodeproject.org/",
            retrieval_date="2026-10-23",
            file_md5="00000000000000000000000000000000",
            peak_chromosome="chr12"
        )
        res_valid = validator.validate(rec_cdk4_valid)
        assert not any("transcript version mismatch" in f for f in res_valid.failures)

        rec_cdk4_invalid = EvidenceRecord(
            gene_symbol="CDK4",
            ncbi_gene_id="1019",
            ensembl_gene_id="ENSG00000135446",
            transcript_id="ENST00000257904",
            transcript_version="99",
            genome_assembly="GRCh38",
            chromosome="chr12",
            strand="-",
            canonical_tss=57752310,
            promoter_upstream=2000,
            promoter_downstream=500,
            promoter_start=57751810,
            promoter_end=57754310,
            encode_file_accession="ENCFF000TEST",
            encode_experiment_accession="ENCSR000TEST",
            experiment_target="CDK4",
            biosample="K562",
            treatment_context="Untreated",
            genetic_modification_context="None",
            source_peak_start=57752000,
            source_peak_end=57752500,
            source_peak_signal=50.0,
            source_peak_summit=57752250,
            overlap_bp=500,
            overlap_percentage=100.0,
            evidence_source_url="https://www.encodeproject.org/",
            retrieval_date="2026-10-23",
            file_md5="00000000000000000000000000000000",
            peak_chromosome="chr12"
        )
        res_invalid = validator.validate(rec_cdk4_invalid)
        assert any(
            "Gene 'CDK4' transcript version mismatch: Authoritative MANE Select version is '11', but record claims '99'."
            in f for f in res_invalid.failures
        )
