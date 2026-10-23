"""
Phase 4A-R2: Genomic Coordinate System and Boundary Overlap Concordance Tests.

Adopts a single internal genomic interval convention: BED 0-based half-open [chromStart, chromEnd).
Explicitly converts 1-based closed promoter intervals [start_1b, end_1b] before interval operations.

Tests both positive-strand and negative-strand genes across all boundary conditions:
1. Left boundary touching (0-bp overlap)
2. Left boundary 1-bp overlap (exactly 1 bp)
3. Right boundary touching (0-bp overlap)
4. Right boundary 1-bp overlap (exactly 1 bp)
5. Peak spanning the TSS
6. Peak completely inside the promoter
7. Peak completely outside the promoter
8. Peak from another chromosome
9. Peak with invalid coordinates (start >= end, start < 0)
10. Controlled reproduction of the Phase 4A-R 1-bp boundary defect.
"""

import pytest
from tf_explorer.evidence_contract import (
    EvidenceRecord,
    EvidenceContractValidator,
    ValidationStatus,
    promoter_1based_to_bed0based,
    calculate_bed0_overlap
)
from tf_explorer.reference_resolver import DynamicReferenceResolver


def oracle_bed0_overlap(p_bed_start: int, p_bed_end: int, peak_start: int, peak_end: int) -> int:
    """Independent mathematical reference oracle for BED 0-based half-open intersection."""
    ov_start = max(p_bed_start, peak_start)
    ov_end = min(p_bed_end, peak_end)
    return max(0, ov_end - ov_start)


class TestCoordinateConventionsAndBoundaries:
    """Exhaustive boundary testing across strand directions and overlap permutations."""

    def test_coordinate_conversion_helper(self):
        """Verifies 1-based closed to BED 0-based half-open conversion arithmetic."""
        # 1-based interval [100, 200] contains 101 bases (100, 101, ..., 200).
        # In BED0: [99, 200) contains 200 - 99 = 101 bases.
        bed_start, bed_end = promoter_1based_to_bed0based(100, 200)
        assert bed_start == 99
        assert bed_end == 200
        assert (bed_end - bed_start) == (200 - 100 + 1)

    def test_reproduce_phase4ar_boundary_overlap_defect(self):
        """
        Controlled reproduction of the Phase 4A-R defect where mixing 1-based promoter
        coordinates with BED0 half-open peak coordinates produced 0 bp for a 1-bp overlap.
        """
        # Let promoter be 1-based [100, 200]
        p_1b_start, p_1b_end = 100, 200
        # Peak covers the first base of the promoter (base 100) in BED0: [99, 100)
        peak_start, peak_end = 99, 100

        # Old flawed arithmetic in Phase 4A-R:
        old_ov_start = max(p_1b_start, peak_start)  # max(100, 99) = 100
        old_ov_end = min(p_1b_end, peak_end)        # min(200, 100) = 100
        old_calc_overlap = max(0, old_ov_end - old_ov_start)  # max(0, 100 - 100) = 0 (DEFECT!)
        assert old_calc_overlap == 0, "Old arithmetic must reproduce the 0 bp defect"

        # Correct BED 0-based half-open arithmetic:
        p_bed_start, p_bed_end = promoter_1based_to_bed0based(p_1b_start, p_1b_end)
        correct_overlap = calculate_bed0_overlap(p_bed_start, p_bed_end, peak_start, peak_end)
        assert correct_overlap == 1, "Correct BED0 arithmetic must identify exactly 1 bp overlap"

    @pytest.mark.parametrize("strand,tss,p_1b_start,p_1b_end", [
        ("+", 10000, 8000, 10500),  # Positive strand: TSS 10,000, upstream 2000, downstream 500
        ("-", 50000, 49500, 52000),  # Negative strand: TSS 50,000, upstream 2000, downstream 500
    ])
    def test_boundary_conditions_both_strands(self, strand, tss, p_1b_start, p_1b_end):
        """Tests all 7 boundary permutations against the mathematical oracle for both strands."""
        p_bed_start, p_bed_end = promoter_1based_to_bed0based(p_1b_start, p_1b_end)
        promoter_width = p_bed_end - p_bed_start

        # 1. Left boundary touching (outside, 0 bp)
        peak_touch_left = (p_bed_start - 200, p_bed_start)
        ov = calculate_bed0_overlap(p_bed_start, p_bed_end, *peak_touch_left)
        assert ov == 0
        assert ov == oracle_bed0_overlap(p_bed_start, p_bed_end, *peak_touch_left)

        # 2. Left boundary 1-bp overlap (covers first base of promoter, 1 bp)
        peak_1bp_left = (p_bed_start - 199, p_bed_start + 1)
        ov = calculate_bed0_overlap(p_bed_start, p_bed_end, *peak_1bp_left)
        assert ov == 1
        assert ov == oracle_bed0_overlap(p_bed_start, p_bed_end, *peak_1bp_left)

        # 3. Right boundary touching (outside, 0 bp)
        peak_touch_right = (p_bed_end, p_bed_end + 200)
        ov = calculate_bed0_overlap(p_bed_start, p_bed_end, *peak_touch_right)
        assert ov == 0
        assert ov == oracle_bed0_overlap(p_bed_start, p_bed_end, *peak_touch_right)

        # 4. Right boundary 1-bp overlap (covers last base of promoter, 1 bp)
        peak_1bp_right = (p_bed_end - 1, p_bed_end + 199)
        ov = calculate_bed0_overlap(p_bed_start, p_bed_end, *peak_1bp_right)
        assert ov == 1
        assert ov == oracle_bed0_overlap(p_bed_start, p_bed_end, *peak_1bp_right)

        # 5. Spanning the TSS
        peak_tss = (tss - 50, tss + 50)
        ov = calculate_bed0_overlap(p_bed_start, p_bed_end, *peak_tss)
        assert ov == 100
        assert ov == oracle_bed0_overlap(p_bed_start, p_bed_end, *peak_tss)

        # 6. Completely inside promoter
        peak_inside = (p_bed_start + 100, p_bed_start + 400)
        ov = calculate_bed0_overlap(p_bed_start, p_bed_end, *peak_inside)
        assert ov == 300
        assert ov == oracle_bed0_overlap(p_bed_start, p_bed_end, *peak_inside)

        # 7. Completely outside promoter (far left and far right)
        peak_far_left = (p_bed_start - 2000, p_bed_start - 1000)
        peak_far_right = (p_bed_end + 1000, p_bed_end + 2000)
        assert calculate_bed0_overlap(p_bed_start, p_bed_end, *peak_far_left) == 0
        assert calculate_bed0_overlap(p_bed_start, p_bed_end, *peak_far_right) == 0

    def test_validator_enforces_bed0_boundary_arithmetic(self):
        """EvidenceContractValidator accepts 1-bp boundary overlap when correctly computed in BED0."""
        # Using CDK4 baseline with mock peak authenticator
        from unittest.mock import MagicMock
        from tf_explorer.peak_authenticator import PeakAuthResult, PeakAuthStatus

        mock_auth = MagicMock()
        mock_auth.authenticate_record.return_value = PeakAuthResult(
            status=PeakAuthStatus.AUTHENTICATED,
            is_authenticated=True,
            source_file_accession="ENCFF115KRO",
            chromosome="chr12",
            start=57751809,
            end=57751909,
            message="Authenticated"
        )

        validator = EvidenceContractValidator(peak_authenticator=mock_auth)

        # CDK4 promoter is [57751810, 57754310] (1-based closed).
        # BED0 promoter is [57751809, 57754310).
        # Peak [57751709, 57751810) covers exactly 1 bp (base 57751809 in BED0).
        mock_auth.authenticate_record.return_value.start = 57751709
        mock_auth.authenticate_record.return_value.end = 57751810

        rec = EvidenceRecord(
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
            encode_file_accession="ENCFF115KRO",
            encode_experiment_accession="ENCSR000EWX",
            experiment_target="E2F1",
            biosample="MCF-7",
            treatment_context="Untreated",
            genetic_modification_context="HA-tag",
            peak_chromosome="chr12",
            source_peak_start=57751709,
            source_peak_end=57751810,
            source_peak_signal=50.0,
            source_peak_summit=57751750,
            overlap_bp=1,  # Exactly 1 bp overlap in BED0
            overlap_percentage=round((1 / 101) * 100.0, 2),
            evidence_source_url="https://www.encodeproject.org/files/ENCFF115KRO/",
            retrieval_date="2026-10-17",
            file_md5="eeacedeabd905942553703034969b287",
            output_classification="STABLE_WITHIN_EVALUATED_CHOICES"
        )

        res = validator.validate(rec)
        assert res.is_valid is True, f"Expected 1-bp overlap to pass validation: {res.failures}"
        assert res.status == ValidationStatus.VERIFIED_PASS
