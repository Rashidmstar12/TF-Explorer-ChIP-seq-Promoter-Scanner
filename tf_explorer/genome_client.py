import requests
import mygene
import logging
from typing import Tuple, Optional

# Configure logging
logger = logging.getLogger(__name__)

ENSEMBL_REST_URL = "https://rest.ensembl.org"

def get_gene_coordinates(gene_name: str, genome: str = "hg38") -> Optional[Tuple[str, int, int, str]]:
    """
    Fetches the TSS coordinates for a given gene using MyGene.info.
    Returns (chrom, start, end, strand) of the gene/transcript.
    
    Note: 'start' and 'end' here are the gene boundaries. 
    For TSS:
      - If strand is +, TSS is start.
      - If strand is -, TSS is end.
    """
    logger.info(f"Fetching coordinates for {gene_name} ({genome})...")
    mg = mygene.MyGeneInfo()
    
    # Map genome build names if necessary. MyGene uses 'hg38' etc.
    # species='human' defaults to hg38 usually.
    
    try:
        # Query for the gene
        # scopes="symbol", fields="genomic_pos"
        res = mg.query(gene_name, scopes="symbol", species="human", fields="genomic_pos,symbol")
        
        if not res or not res.get('hits') or len(res['hits']) == 0:
            logger.error(f"Gene {gene_name} not found.")
            return None
            
        # Iterate through hits to find exact symbol match
        target_hit = None
        for hit in res['hits']:
            if hit.get('symbol') == gene_name:
                target_hit = hit
                break
        
        # If no exact match, use the first one (fallback)
        if not target_hit:
            logger.warning(f"No exact match for symbol '{gene_name}'. Using first hit: {res['hits'][0].get('symbol')}")
            target_hit = res['hits'][0]
            
        # genomic_pos can be a list (multiple locations) or dict
        gpos = target_hit.get('genomic_pos')
        
        if isinstance(gpos, list):
            # If multiple, try to find the one matching the requested genome (though mygene usually returns current assembly)
            # For simplicity, take the first one or the one that looks like a primary assembly.
            # Usually mygene returns hg38 for human.
            gpos = gpos[0]
            
        if not gpos:
            logger.error(f"No genomic position found for {gene_name}.")
            return None
            
        chrom = str(gpos['chr'])
        start = int(gpos['start'])
        end = int(gpos['end'])
        strand = str(gpos['strand']) # 1 or -1
        
        # Convert strand to +/-
        strand_symbol = "+" if str(strand) == "1" else "-"
        
        logger.info(f"Found {gene_name} at chr{chrom}:{start}-{end} ({strand_symbol})")
        return chrom, start, end, strand_symbol
        
    except Exception as e:
        logger.error(f"Error fetching gene coordinates: {e}")
        return None

def get_promoter_sequence(chrom: str, tss: int, strand: str, up_bp: int, down_bp: int, genome: str = "hg38") -> Optional[str]:
    """
    Fetches the promoter sequence from Ensembl REST API.
    
    chrom: Chromosome name (e.g., "1", "X") - Ensembl expects just the number/letter, no 'chr' prefix usually, 
           but we should handle both.
    tss: Transcription Start Site coordinate.
    strand: "+" or "-"
    up_bp: Base pairs upstream of TSS.
    down_bp: Base pairs downstream of TSS.
    """
    # Calculate region
    if strand == "+":
        # Upstream is < TSS, Downstream is > TSS
        # Region: [TSS - up, TSS + down]
        seq_start = tss - up_bp
        seq_end = tss + down_bp
    else:
        # Strand is -
        # Upstream is > TSS, Downstream is < TSS (genomic coordinates)
        # But usually "upstream" means 5' of the gene.
        # For - strand gene, 5' end is at the larger genomic coordinate (end of gene).
        # So "upstream" is > TSS.
        # Region: [TSS - down, TSS + up]
        # Wait, let's be precise.
        # If gene is on -, TSS is at genomic 'end'.
        # 2000bp upstream means genomic coordinates: TSS to TSS + 2000.
        # 500bp downstream means genomic coordinates: TSS - 500 to TSS.
        # So region is [TSS - down, TSS + up]
        seq_start = tss - down_bp
        seq_end = tss + up_bp

    # Ensembl API expects coords without 'chr' prefix usually
    clean_chrom = chrom.replace("chr", "")
    
    logger.info(f"Fetching sequence for chr{clean_chrom}:{seq_start}-{seq_end}...")
    
    ext = f"/sequence/region/human/{clean_chrom}:{seq_start}..{seq_end}:1" # :1 for + strand of reference
    # We want the sequence on the coding strand of the gene?
    # Usually promoter analysis is done on the coding strand (relative to gene).
    # If strand is -, we should request strand=-1 from Ensembl to get the reverse complement (coding sequence).
    
    api_strand = 1 if strand == "+" else -1
    ext = f"/sequence/region/human/{clean_chrom}:{seq_start}..{seq_end}:{api_strand}"
    
    headers = {"Content-Type": "text/plain"}
    
    try:
        r = requests.get(ENSEMBL_REST_URL + ext, headers=headers)
        
        if not r.ok:
            logger.error(f"Ensembl API failed: {r.text}")
            r.raise_for_status()
            
        return r.text
    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching sequence: {e}")
        return None
