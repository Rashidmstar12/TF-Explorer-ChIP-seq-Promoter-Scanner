import pytest
import numpy as np
from Bio import motifs
from Bio.Seq import Seq
import tf_explorer.motifs as tem

def test_motif_strand_detection_and_no_duplicates():
    """
    Test that scan_promoter_with_jaspar:
    1. Accurately detects forward strand (+) matches.
    2. Accurately detects reverse strand (-) matches.
    3. Does not double count matches.
    """
    # Create synthetic non-palindromic motif: ACCC (length 4, reverse complement GGGT)
    instances = [Seq("ACCC")]
    m = motifs.create(instances)
    m.name = "TEST_TF"
    
    # Test sequence: TTTACCCAATTTGGGTAA (length 18)
    # Forward match: ACCC at 3..7
    # Reverse match: GGGT at 12..16 (since RC of GGGT is ACCC)
    test_seq = "TTTACCCAATTTGGGTAA"
    
    original_get = tem.get_jaspar_matrix
    tem.get_jaspar_matrix = lambda mid: m
    try:
        hits = tem.scan_promoter_with_jaspar(test_seq, ["TEST_ID"], threshold=2.0)
        
        # Must have exactly 2 hits: one on +, one on -
        assert len(hits) == 2, f"Expected exactly 2 hits, got {len(hits)}"
        
        hit_fwd = [h for h in hits if h["strand"] == "+"]
        hit_rev = [h for h in hits if h["strand"] == "-"]
        
        assert len(hit_fwd) == 1, f"Expected 1 forward hit, got {len(hit_fwd)}"
        assert len(hit_rev) == 1, f"Expected 1 reverse hit, got {len(hit_rev)}"
        
        # Verify forward hit coordinates and sequence
        assert hit_fwd[0]["start"] == 3
        assert hit_fwd[0]["end"] == 7
        assert hit_fwd[0]["sequence"] == "ACCC"
        
        # Verify reverse hit coordinates and sequence
        assert hit_rev[0]["start"] == 12
        assert hit_rev[0]["end"] == 16
        assert hit_rev[0]["sequence"] == "ACCC"
        
    finally:
        tem.get_jaspar_matrix = original_get

def test_palindromic_motif_both_strands():
    """
    Palindromic motif (e.g. CGCG) matches both orientations at the same locus.
    Verifies exactly one hit per strand orientation.
    """
    pal_m = motifs.create([Seq("CGCG")])
    pal_m.name = "PAL_TF"
    
    test_seq = "TTTCGCGTTT" # length 10, match at 3..7
    
    original_get = tem.get_jaspar_matrix
    tem.get_jaspar_matrix = lambda mid: pal_m
    try:
        hits = tem.scan_promoter_with_jaspar(test_seq, ["PAL_ID"], threshold=2.0)
        assert len(hits) == 2, f"Expected 2 hits (1 per strand) for palindrome, got {len(hits)}"
        strands = {h["strand"] for h in hits}
        assert strands == {"+", "-"}
        for h in hits:
            assert h["start"] == 3
            assert h["end"] == 7
    finally:
        tem.get_jaspar_matrix = original_get
