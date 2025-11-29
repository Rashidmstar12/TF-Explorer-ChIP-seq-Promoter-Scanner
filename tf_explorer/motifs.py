import requests
import logging
from typing import List, Dict, Tuple
from Bio import motifs
from Bio.Seq import Seq
from io import StringIO

# Configure logging
logger = logging.getLogger(__name__)

JASPAR_API_URL = "https://jaspar.genereg.net/api/v1/matrix"

# Common TFs to JASPAR IDs mapping (for convenience/defaults)
COMMON_TFS = {
    "CTCF": "MA0139.1",
    "SP1": "MA0079.3",
    "E2F1": "MA0024.3",
    "YY1": "MA0095.2",
    "TP53": "MA0106.3",
    "MYC": "MA0147.3",
    "CREB1": "MA0018.3",
    "CREB": "MA0018.3"
}

def get_jaspar_matrix(matrix_id: str):
    """
    Fetches a PFM from JASPAR API and returns a Bio.motifs object.
    """
    logger.info(f"Fetching JASPAR matrix {matrix_id}...")
    try:
        # Fetch PFM format
        url = f"{JASPAR_API_URL}/{matrix_id}.pfm"
        response = requests.get(url)
        response.raise_for_status()
        
        # Parse with Bio.motifs
        # Bio.motifs.read expects a handle
        handle = StringIO(response.text)
        m = motifs.read(handle, "jaspar")
        return m
    except Exception as e:
        logger.error(f"Failed to fetch/parse matrix {matrix_id}: {e}")
        return None

def scan_promoter_with_jaspar(promoter_seq: str, jaspar_ids: List[str], threshold: float = 8.0) -> List[Dict]:
    """
    Scans the promoter sequence with the given JASPAR matrices.
    Returns a list of hits.
    """
    if jaspar_ids is None:
        return []
    if len(jaspar_ids) == 0:
        return []
        
    logger.info(f"Scanning promoter ({len(promoter_seq)} bp) with {len(jaspar_ids)} matrices...")
    
    hits = []
    seq_obj = Seq(promoter_seq)
    
    for mid in jaspar_ids:
        m = get_jaspar_matrix(mid)
        if not m:
            continue
            
        # Convert to PWM (PSSM in Biopython terms)
        # Add pseudocounts to avoid log(0)
        pwm = m.counts.normalize(pseudocounts=0.5).log_odds()
        
        # Search
        # pssm.search returns (position, score) tuples
        # It searches both strands by default? No, usually just forward.
        # We should check both strands or just the provided one.
        # Promoter analysis usually considers both strands as TFs can bind either way.
        
        # Search forward
        for pos, score in pwm.search(seq_obj, threshold=threshold):
            hits.append({
                "tf_name": m.name,
                "jaspar_id": mid,
                "strand": "+",
                "start": pos,
                "end": pos + len(m),
                "score": score,
                "sequence": str(seq_obj[pos:pos+len(m)])
            })
            
        # Search reverse
        rc_seq = seq_obj.reverse_complement()
        for pos, score in pwm.search(rc_seq, threshold=threshold):
            # Position on RC needs to be mapped back to forward if we want genomic coords relative to TSS?
            # Or just report it as relative to 5' of the provided sequence.
            # Let's report relative to the start of the provided sequence.
            # If match is at pos on RC, it means it's at len(seq) - pos - len(m) on forward.
            
            fwd_start = len(promoter_seq) - pos - len(m)
            fwd_end = len(promoter_seq) - pos
            
            hits.append({
                "tf_name": m.name,
                "jaspar_id": mid,
                "strand": "-",
                "start": fwd_start,
                "end": fwd_end,
                "score": score,
                "sequence": str(rc_seq[pos:pos+len(m)]) # This is the sequence on the - strand
            })
            
    return hits
