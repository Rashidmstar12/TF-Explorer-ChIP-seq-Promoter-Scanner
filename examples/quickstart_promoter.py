"""
Quickstart Example: Promoter Boundary & Coordinate Calculation
Demonstrates strand-aware promoter coordinate expansion for human genes.
"""

from tf_explorer.analysis import calculate_promoter_coordinates

def main():
    print("=== TF-Explorer: Strand-Aware Promoter Coordinate Expansion ===")
    
    # Case 1: Positive-strand gene (e.g., CTCF on chr16, canonical TSS: 67,562,526)
    gene_pos = "CTCF"
    tss_pos = 67562526
    strand_pos = "+"
    p_start, p_end = calculate_promoter_coordinates(tss_pos, strand_pos, upstream=2000, downstream=500)
    print(f"Gene: {gene_pos} ({strand_pos} strand, TSS: {tss_pos})")
    print(f"  Promoter interval (-2000 bp / +500 bp): chr16:{p_start}-{p_end} (Length: {p_end - p_start} bp)\n")
    
    # Case 2: Negative-strand gene (e.g., CDK4 on chr12, canonical TSS: 57,752,310)
    gene_neg = "CDK4"
    tss_neg = 57752310
    strand_neg = "-"
    n_start, n_end = calculate_promoter_coordinates(tss_neg, strand_neg, upstream=2000, downstream=500)
    print(f"Gene: {gene_neg} ({strand_neg} strand, TSS: {tss_neg})")
    print(f"  Promoter interval (-2000 bp / +500 bp): chr12:{n_start}-{n_end} (Length: {n_end - n_start} bp)")
    print(f"  Note: On negative strand, upstream expands to higher genomic coordinate.")

if __name__ == "__main__":
    main()
