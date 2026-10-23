import pytest
import gzip
import tempfile
import os
import pandas as pd
from tf_explorer.analysis import find_overlaps_for_experiment

def test_strict_overlap_and_nonoverlapping_separation():
    """
    Verifies that:
    1. Peaks strictly overlapping the promoter are returned in hits_strict with overlap=True.
    2. Peaks outside the promoter are NOT in hits_strict.
    3. The closest non-overlapping peak is identified with overlap=False.
    4. Zero strict overlaps yields hits_strict == [].
    """
    # Create synthetic BED file with two peaks:
    # 1. Distant peak: chr1:50000-50500 (distance = 40,250 bp)
    # 2. Strictly overlapping peak: chr1:9500-9800 (inside promoter [8000, 10500))
    with tempfile.NamedTemporaryFile(suffix=".bed.gz", delete=False) as tmp:
        bed_path = tmp.name
        
    try:
        with gzip.open(bed_path, "wt") as f:
            # col0: chrom, col1: start, col2: end, col3: name, col4: score, col5: strand, col6: signalValue
            f.write("chr1\t50000\t50500\tPEAK_DISTANT\t100\t.\t10.5\n")
            f.write("chr1\t9500\t9800\tPEAK_PROMOTER\t200\t.\t25.0\n")
            
        exp_meta = {'dataset_accession': 'ENCSR000TEST', 'file_accession': 'ENCFF000TEST', 'biosample': 'K562'}
        
        # Test 1: Window [8000, 10500) overlaps PEAK_PROMOTER
        hits_strict, hits_loose, total, closest = find_overlaps_for_experiment(
            bed_path, chrom="1",
            p_start_strict=8000, p_end_strict=10500,
            p_start_loose=5000, p_end_loose=15000,
            tss=10000, strand="+", exp=exp_meta, tf="TEST_TF"
        )
        
        assert len(hits_strict) == 1, f"Expected 1 strict hit, got {len(hits_strict)}"
        assert hits_strict[0]["peak_start"] == 9500
        assert hits_strict[0]["overlap"] is True
        
        # Test 2: Window [8000, 9000) does NOT overlap either peak
        hits_strict_2, hits_loose_2, total_2, closest_2 = find_overlaps_for_experiment(
            bed_path, chrom="1",
            p_start_strict=8000, p_end_strict=9000,
            p_start_loose=7000, p_end_loose=9200,
            tss=8500, strand="+", exp=exp_meta, tf="TEST_TF"
        )
        
        assert len(hits_strict_2) == 0, f"Expected 0 strict hits, got {len(hits_strict_2)}"
        assert closest_2 is not None, "Closest peak should be identified"
        assert closest_2["peak_start"] == 9500
        assert closest_2["overlap"] is False
        
        # Verify that downstream pipeline logic separates them:
        encode_hits = []
        nearest_nonoverlapping_peaks = []
        if hits_strict_2:
            encode_hits.extend(hits_strict_2)
        elif closest_2:
            nearest_nonoverlapping_peaks.append(closest_2)
            
        assert len(encode_hits) == 0, "Strict encode_hits must be empty when 0 peaks overlap promoter"
        assert len(nearest_nonoverlapping_peaks) == 1, "Nearest non-overlapping peak must be stored in separate list"
        
        df_encode = pd.DataFrame(encode_hits)
        assert len(df_encode) == 0, "Summary peak count must be 0"
        
    finally:
        if os.path.exists(bed_path):
            os.remove(bed_path)
