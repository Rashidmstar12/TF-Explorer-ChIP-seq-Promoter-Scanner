import requests
import logging
from typing import List, Dict, Tuple
from Bio import motifs
from Bio.Seq import Seq
from io import StringIO

# Configure logging
logger = logging.getLogger(__name__)

JASPAR_API_URL = "https://jaspar.elixir.no/api/v1/matrix"
JASPAR_FALLBACK_URL = "https://jaspar.genereg.net/api/v1/matrix"
JASPAR_CACHE_DIR = "jaspar_cache"

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
    Fetches a PFM from JASPAR API (with local file caching) and returns a Bio.motifs object.
    """
    logger.info(f"Fetching JASPAR matrix {matrix_id}...")
    import os
    os.makedirs(JASPAR_CACHE_DIR, exist_ok=True)
    cache_path = os.path.join(JASPAR_CACHE_DIR, f"{matrix_id}.pfm")

    # 1. Check local cache
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as f:
                content = f.read()
                handle = StringIO(content)
                m = motifs.read(handle, "jaspar")
                if m is not None and m.length > 0:
                    return m
        except Exception as e_cache:
            logger.warning(f"Purging corrupted cached matrix {matrix_id}: {e_cache}")
            try:
                os.remove(cache_path)
            except OSError:
                pass

    # 2. Fetch from active API
    for api_base in [JASPAR_API_URL, JASPAR_FALLBACK_URL]:
        try:
            url = f"{api_base}/{matrix_id}.pfm"
            response = requests.get(url, timeout=10)
            if response.ok and response.text.strip().startswith(">"):
                handle = StringIO(response.text)
                m = motifs.read(handle, "jaspar")
                if m is not None and m.length > 0:
                    # Validated; save to cache atomically
                    with open(cache_path, "w", encoding="utf-8") as f:
                        f.write(response.text)
                    return m
        except Exception as e:
            logger.warning(f"Failed fetching/parsing {matrix_id} from {api_base}: {e}")

    logger.error(f"Failed to fetch/parse matrix {matrix_id} from all sources.")
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
    seq_obj = Seq(promoter_seq.upper())
    seq_len = len(promoter_seq)
    
    for mid in jaspar_ids:
        m = get_jaspar_matrix(mid)
        if not m:
            continue
            
        # Convert to PWM (PSSM in Biopython terms)
        # Add pseudocounts to avoid log(0)
        pwm = m.counts.normalize(pseudocounts=0.5).log_odds()
        motif_len = len(m)
        
        # Biopython pssm.search evaluates both strands natively when both=True (default).
        # Positive positions (pos >= 0) indicate forward strand matches at [pos, pos + motif_len).
        # Negative positions (pos < 0) indicate reverse strand matches where pos = neg_pos - seq_len.
        # The 0-based sequence offset for reverse strand matches is therefore: seq_len + pos.
        for pos, score in pwm.search(seq_obj, threshold=threshold):
            if pos >= 0:
                start = int(pos)
                end = start + motif_len
                strand = "+"
                hit_seq = str(seq_obj[start:end])
            else:
                start = int(seq_len + pos)
                end = start + motif_len
                strand = "-"
                hit_seq = str(seq_obj[start:end].reverse_complement())
                
            if 0 <= start <= end <= seq_len:
                hits.append({
                    "tf_name": getattr(m, "name", mid) or mid,
                    "jaspar_id": mid,
                    "strand": strand,
                    "start": start,
                    "end": end,
                    "score": float(score),
                    "sequence": hit_seq
                })
                
    hits.sort(key=lambda x: (x["start"], x["end"], x["strand"]))
    return hits
