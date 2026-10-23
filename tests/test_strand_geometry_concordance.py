"""
Authoritative Strand and Genomic Coordinate Geometry Concordance Tests.

Validates that:
1. Negative-strand transcripts (CDK4, TFF1) use authoritative (-) promoter arithmetic.
2. Positive-strand transcripts (MYC) use authoritative (+) promoter arithmetic.
3. Strand inversion is detected and rejected (tests fail if strand logic is reversed).
4. Authentic peaks bracket true biological TSS coordinates without synthetic distortion.
"""

import pytest
from tf_explorer.robustness import RobustnessAnalyzer, RobustnessClassification


class TestStrandGeometryConcordance:
    @pytest.fixture
    def analyzer(self):
        return RobustnessAnalyzer()

    def test_cdk4_negative_strand_promoter_geometry(self, analyzer):
        """
        CDK4 (NCBI:1019, Ensembl:ENSG00000135446) is on chr12 NEGATIVE strand.
        Canonical MANE Select TSS: 57,752,310 (ENST00000257904.11).
        Promoter window (-2000/+500) evaluates to [57751810, 57754310].
        """
        tss = 57752310
        strand = "-"
        p_start, p_end = analyzer.calculate_promoter_coordinates(tss, strand, 2000, 500)
        assert p_start == 57751810
        assert p_end == 57754310
        assert p_end - p_start == 2500

        # Overlap with authentic MCF-7 E2F1 peak (ENCFF115KRO: [57752169, 57752512])
        is_ov, ov_bp = analyzer.check_peak_overlap(57752169, 57752512, p_start, p_end)
        assert is_ov is True
        assert ov_bp == (57752512 - 57752169)  # 343 bp full overlap

        # Inversion check: positive strand arithmetic produces displaced coordinates
        p_wrong_start, p_wrong_end = analyzer.calculate_promoter_coordinates(tss, "+", 2000, 500)
        assert p_wrong_start == 57750310
        assert p_wrong_end == 57752810
        assert p_start != p_wrong_start

    def test_tff1_authoritative_negative_strand_geometry(self, analyzer):
        """
        TFF1 (NCBI Gene ID: 7031, Ensembl: ENSG00000160182) is on chr21 NEGATIVE strand.
        Canonical MANE Select TSS: 42,366,535.
        Promoter window (-2000/+500) must evaluate to [42366035, 42368535].
        """
        tss = 42366535
        strand = "-"
        p_start, p_end = analyzer.calculate_promoter_coordinates(tss, strand, 2000, 500)
        assert p_start == 42366035, f"Expected 42366035, got {p_start}"
        assert p_end == 42368535, f"Expected 42368535, got {p_end}"
        assert p_end - p_start == 2500

        # Overlap with authentic ESR1 peak (ENCFF185UUT: [42366395, 42367086])
        is_ov, ov_bp = analyzer.check_peak_overlap(42366395, 42367086, p_start, p_end)
        assert is_ov is True
        assert ov_bp == (42367086 - 42366395)  # 691 bp full overlap

        # Anti-regression check: If mistakenly treated as positive strand (+)
        p_pos_start, p_pos_end = analyzer.calculate_promoter_coordinates(tss, "+", 2000, 500)
        assert p_pos_start == 42364535
        assert p_pos_end == 42367035
        # The positive-strand window truncates the upstream flank of the peak
        assert p_start != p_pos_start

    def test_myc_authoritative_positive_strand_geometry(self, analyzer):
        """
        MYC (Ensembl: ENSG00000136997) is on chr8 POSITIVE strand.
        Canonical MANE Select TSS: 127,736,231 (ENST00000621592.8).
        Promoter window (-2000/+500) must evaluate to [127734231, 127736731].
        """
        tss = 127736231
        strand = "+"
        p_start, p_end = analyzer.calculate_promoter_coordinates(tss, strand, 2000, 500)
        assert p_start == 127734231, f"Expected 127734231, got {p_start}"
        assert p_end == 127736731, f"Expected 127736731, got {p_end}"
        assert p_end - p_start == 2500

        # Overlap with authentic E2F1 peak in HeLa-S3 (ENCFF004ZKN: [127735986, 127736362])
        is_ov, ov_bp = analyzer.check_peak_overlap(127735986, 127736362, p_start, p_end)
        assert is_ov is True
        assert ov_bp == (127736362 - 127735986)  # 376 bp full overlap

        # Inversion check: Negative strand arithmetic on MYC produces wrong coordinates
        p_neg_start, p_neg_end = analyzer.calculate_promoter_coordinates(tss, "-", 2000, 500)
        assert p_neg_start == 127735731
        assert p_neg_end == 127738231
        assert p_start != p_neg_start

    def test_tff1_multiverse_evaluation_negative_strand(self, analyzer):
        """End-to-end evaluation of TFF1 ESR1 in MCF-7 under negative strand geometry."""
        transcripts = [
            {"id": "ENST00000291527.3", "name": "TFF1-201 (Canonical)", "tss": 42366535, "strand": "-", "is_canonical": True}
        ]
        peaks = [{"chrom": "chr21", "start": 42366395, "end": 42367086, "signal_value": 2601.5}]
        exp_sets = [{"experiment_id": "ENCSR463GOT", "quality_tier": "optimal", "peaks": peaks}]

        rep = analyzer.evaluate(
            gene_symbol="TFF1",
            target_tf="ESR1",
            biosample="MCF-7",
            transcripts=transcripts,
            experiment_peak_sets=exp_sets,
            promoter_windows=[(500, 100), (1000, 200), (2000, 500)]
        )
        assert rep.classification == RobustnessClassification.STABLE_WITHIN_EVALUATED_CHOICES
        assert rep.support_ratio == 1.0
        assert rep.supporting_configurations == 3
        assert rep.total_valid_configurations == 3
        assert rep.transcript_support["ENST00000291527.3"] is True
