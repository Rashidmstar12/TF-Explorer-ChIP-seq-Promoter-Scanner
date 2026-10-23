import pytest

def get_promoter_bounds(tss: int, strand: str, up_bp: int, down_bp: int):
    """
    Standard formula implemented in tf_explorer.analysis:
    - Positive strand (+): [tss - up_bp, tss + down_bp)
    - Negative strand (-): [tss - down_bp, tss + up_bp)
    """
    if strand == "+":
        return tss - up_bp, tss + down_bp
    else:
        return tss - down_bp, tss + up_bp

def test_positive_strand_promoter_interval():
    tss = 10000
    up = 2000
    down = 500
    start, end = get_promoter_bounds(tss, "+", up, down)
    assert start == 8000, f"Expected start=8000, got {start}"
    assert end == 10500, f"Expected end=10500, got {end}"
    assert end - start == 2500, f"Expected length=2500, got {end - start}"

def test_negative_strand_promoter_interval():
    tss = 25000
    up = 2000
    down = 500
    start, end = get_promoter_bounds(tss, "-", up, down)
    assert start == 24500, f"Expected start=24500, got {start}"
    assert end == 27000, f"Expected end=27000, got {end}"
    assert end - start == 2500, f"Expected length=2500, got {end - start}"

def test_interval_overlap_half_open_logic():
    """
    Verifies standard half-open interval overlap: max(s1, s2) < min(e1, e2)
    """
    p_start, p_end = 1000, 2000
    
    def overlaps(s, e):
        return max(p_start, s) < min(p_end, e)
    
    # Internal peak
    assert overlaps(1200, 1500) is True
    # Overlapping 5' boundary
    assert overlaps(800, 1200) is True
    # Overlapping 3' boundary
    assert overlaps(1800, 2200) is True
    # Touching 5' boundary [500, 1000): 0bp overlap -> False
    assert overlaps(500, 1000) is False
    # Touching 3' boundary [2000, 2500): 0bp overlap -> False
    assert overlaps(2000, 2500) is False
    # Far upstream [100, 300) -> False
    assert overlaps(100, 300) is False
    # Far downstream [3000, 3500) -> False
    assert overlaps(3000, 3500) is False
    # Zero-length feature at 1500 [1500, 1500) -> False
    assert overlaps(1500, 1500) is False
