import argparse
import sys
import os
import random
import logging
from typing import List

from . import analysis, motifs

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_jaspar_ids(input_str: str) -> List[str]:
    """
    Parses a comma-separated string of JASPAR IDs or TF names.
    Resolves names to IDs using motifs.COMMON_TFS.
    """
    if not input_str:
        return []
        
    ids = []
    for item in input_str.split(','):
        item = item.strip()
        if item in motifs.COMMON_TFS:
            ids.append(motifs.COMMON_TFS[item])
        else:
            ids.append(item)
    return ids

def main():
    parser = argparse.ArgumentParser(description="tf-explorer: Analyze TF binding sites for ANY gene.")
    
    # Input
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene", type=str, help="Target gene symbol (e.g., TP53)")
    group.add_argument("--batch-genes", type=str, help="Path to a text file containing a list of genes (one per line)")
    
    parser.add_argument("--tf-list", type=str, default="", help="Comma-separated list of TFs to search in ENCODE (e.g., E2F1,MYC)")
    parser.add_argument("--jaspar-ids", type=str, default="", help="Comma-separated list of JASPAR Matrix IDs or TF names (e.g., MA0139.1,CTCF)")
    
    # Parameters
    parser.add_argument("--genome", type=str, default="hg38", help="Genome assembly (default: hg38)")
    parser.add_argument("--promoter-up", type=int, default=2000, help="Base pairs upstream of TSS (default: 2000)")
    parser.add_argument("--promoter-down", type=int, default=500, help="Base pairs downstream of TSS (default: 500)")
    parser.add_argument("--threshold", type=float, default=8.0, help="PWM score threshold (default: 8.0)")
    
    # Output & Features
    parser.add_argument("--out", type=str, default="results", help="Output directory")
    parser.add_argument("--bed-output", action="store_true", help="Generate BED file of overlaps")
    parser.add_argument("--plot-track", action="store_true", help="Generate promoter track plot")
    parser.add_argument("--random-seed", type=int, default=42, help="Random seed for reproducibility")
    
    args = parser.parse_args()
    
    # Set seed
    random.seed(args.random_seed)
    
    # Parse lists
    tf_list = [x.strip() for x in args.tf_list.split(',')] if args.tf_list else []
    
    # Resolve JASPAR IDs
    jaspar_ids = parse_jaspar_ids(args.jaspar_ids)
    
    # Default JASPAR IDs if none provided (User request: "INCLUDE 2 OR 3 TF PREDICTION ALSO")
    if not jaspar_ids:
        logger.info("No JASPAR IDs provided. Using default set (CTCF, SP1, TP53).")
        jaspar_ids = [motifs.COMMON_TFS["CTCF"], motifs.COMMON_TFS["SP1"], motifs.COMMON_TFS["TP53"]]
        
    # Genes to process
    genes = []
    if args.gene:
        genes.append(args.gene)
    elif args.batch_genes:
        if os.path.exists(args.batch_genes):
            with open(args.batch_genes, 'r') as f:
                genes = [line.strip() for line in f if line.strip()]
        else:
            logger.error(f"Batch file {args.batch_genes} not found.")
            sys.exit(1)
            
    logger.info(f"Processing {len(genes)} gene(s)...")
    
    for gene in genes:
        # Create gene-specific output subdir if batch processing?
        # Or just put everything in --out?
        # If batch, probably better to have subdirs or prefixed files.
        # analysis.py uses gene_name prefix, so flat dir is okay, but subdirs are cleaner.
        
        gene_out_dir = os.path.join(args.out, gene) if len(genes) > 1 else args.out
        
        try:
            analysis.analyze_gene(
                gene_name=gene,
                tf_list=tf_list,
                jaspar_ids=jaspar_ids,
                genome=args.genome,
                promoter_up=args.promoter_up,
                promoter_down=args.promoter_down,
                pwm_threshold=args.threshold,
                output_dir=gene_out_dir,
                bed_output=args.bed_output,
                plot_track=args.plot_track,
                random_seed=args.random_seed
            )
        except Exception as e:
            logger.error(f"Failed to analyze {gene}: {e}")
            continue

if __name__ == "__main__":
    main()
