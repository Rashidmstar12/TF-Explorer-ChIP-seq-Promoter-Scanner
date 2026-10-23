import pytest
from Bio import motifs
from Bio.Seq import Seq
import tf_explorer.motifs as tem

def test_motif_coordinate_bounds_and_intervals():
    """
    Verifies that:
    1. start >= 0 always (no negative local offsets).
    2. end <= len(sequence) always (no out-of-bounds coordinates).
    3. end == start + motif_length.
    4. Boundary matches at index 0 and index len(seq)-len(m) are exact.
    """
    instances = [Seq("ACCC")]
    m = motifs.create(instances)
    m.name = "BOUNDARY_TF"
    
    # Boundary sequence: ACCC at 0..4 (+), and GGGT at 11..15 (-)
    test_seq = "ACCCAAAAAAAGGGTA" # length 16
    
    original_get = tem.get_jaspar_matrix
    tem.get_jaspar_matrix = lambda mid: m
    try:
        hits = tem.scan_promoter_with_jaspar(test_seq, ["BND_ID"], threshold=2.0)
        
        assert len(hits) == 2, f"Expected 2 boundary hits, got {len(hits)}"
        
        for h in hits:
            start = h["start"]
            end = h["end"]
            assert start >= 0, f"Negative start coordinate detected: {start}"
            assert end <= len(test_seq), f"End coordinate {end} exceeds sequence length {len(test_seq)}"
            assert end == start + 4, f"Interval length mismatch: {end - start} != 4"
            
        # First hit starts exactly at 0
        assert hits[0]["start"] == 0
        assert hits[0]["end"] == 4
        assert hits[0]["strand"] == "+"
        
        # Second hit starts at 11, ends at 15
        assert hits[1]["start"] == 11
        assert hits[1]["end"] == 15
        assert hits[1]["strand"] == "-"
        
    finally:
        tem.get_jaspar_matrix = original_get

def test_motif_ambiguous_nucleotides_handled():
    """
    Verifies that ambiguous nucleotides (e.g. N) do not cause crashes or spurious hits.
    """
    instances = [Seq("ACCC")]
    m = motifs.create(instances)
    m.name = "AMB_TF"
    
    test_seq = "NNNACNCAATTTGGNTAA"
    
    original_get = tem.get_jaspar_matrix
    tem.get_jaspar_matrix = lambda mid: m
    try:
        hits = tem.scan_promoter_with_jaspar(test_seq, ["AMB_ID"], threshold=2.0)
        assert len(hits) == 0, f"Expected 0 hits for degenerate sequence, got {len(hits)}"
    finally:
        tem.get_jaspar_matrix = original_get
