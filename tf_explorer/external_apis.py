"""
Pure-Python external API integration module for TF-Explorer.
Provides native, secure HTTP clients for:
- UCSC Genome Browser REST API (phyloP & phastCons conservation scores)
- GTEx Portal API v2 (tissue expression data)
- STRING Database REST API v12 (protein-protein interaction network & partners)

Eliminates all external shell/subprocess dependencies and machine-specific paths.
"""

import io
import logging
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd
import requests

logger = logging.getLogger(__name__)

# Constants and Base URLs
UCSC_API_URL = "https://api.genome.ucsc.edu/getData/track"
GTEx_BASE_URL = "https://gtexportal.org/api/v2"
STRING_BASE_URL = "https://version-12-0.string-db.org/api"
DEFAULT_TIMEOUT = 10  # seconds


# ============================================================================
# UCSC Genome Browser API (Conservation)
# ============================================================================

def get_ucsc_conservation(chrom: str, start: int, end: int, genome: str = "hg38") -> Dict[str, Any]:
    """
    Fetches base-by-base phyloP and phastCons conservation scores from UCSC REST API.
    
    Returns a dictionary formatted for analysis.py:
    {
        "<chrom>:<start>-<end>": {
            "tracks": {
                "phyloP100way": [...],
                "phastCons100way": [...]
            }
        }
    }
    """
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
        
    genome_norm = genome.lower()
    if genome_norm not in ("hg38", "hg19"):
        genome_norm = "hg38"
        
    tracks_to_fetch = ["phyloP100way", "phastCons100way"]
    tracks_data: Dict[str, List[Dict[str, Any]]] = {}
    
    for track in tracks_to_fetch:
        params = {
            "genome": genome_norm,
            "track": track,
            "chrom": chrom,
            "start": start,
            "end": end
        }
        try:
            resp = requests.get(UCSC_API_URL, params=params, timeout=DEFAULT_TIMEOUT)
            if resp.status_code == 200:
                data = resp.json()
                items = data.get(chrom, [])
                if not items and track in data:
                    items = data.get(track, [])
                tracks_data[track] = items if isinstance(items, list) else []
            else:
                logger.warning(f"UCSC API returned status {resp.status_code} for track {track}")
                tracks_data[track] = []
        except Exception as e:
            logger.warning(f"Failed to fetch UCSC conservation ({track}): {e}")
            tracks_data[track] = []
            
    coord_key = f"{chrom}:{start}-{end}"
    return {coord_key: {"tracks": tracks_data}}


# ============================================================================
# GTEx Portal API (Tissue Expression)
# ============================================================================

def resolve_gtex_gencode_id(gene_symbol: str, gencode_version: str = "v39") -> Optional[str]:
    """
    Resolves a standard gene symbol to its versioned GENCODE ID via GTEx reference API.
    """
    url = f"{GTEx_BASE_URL}/reference/gene"
    params = {
        "geneId": gene_symbol.strip(),
        "gencodeVersion": gencode_version,
        "pageSize": 10
    }
    try:
        resp = requests.get(url, params=params, timeout=DEFAULT_TIMEOUT)
        if resp.status_code != 200:
            return None
        data = resp.json()
        items = data.get("data", []) if isinstance(data, dict) else data
        if not items:
            return None
            
        # Match exact symbol case-insensitively
        best = items[0]
        for item in items:
            if item.get("geneSymbol", "").lower() == gene_symbol.strip().lower():
                best = item
                break
        return best.get("gencodeId")
    except Exception as e:
        logger.warning(f"GTEx gene ID resolution failed for {gene_symbol}: {e}")
        return None


def get_gtex_median_expression(gencode_id: str, dataset_id: str = "gtex_v10") -> List[Dict[str, Any]]:
    """
    Fetches median gene expression (TPM) across all tissues for a given GENCODE ID.
    """
    url = f"{GTEx_BASE_URL}/expression/medianGeneExpression"
    params = {
        "gencodeId": gencode_id,
        "datasetId": dataset_id,
        "pageSize": 250
    }
    try:
        resp = requests.get(url, params=params, timeout=DEFAULT_TIMEOUT)
        if resp.status_code != 200:
            return []
        data = resp.json()
        items = data.get("data", []) if isinstance(data, dict) else data
        return items if isinstance(items, list) else []
    except Exception as e:
        logger.warning(f"GTEx median expression query failed for {gencode_id}: {e}")
        return []


def fetch_gtex_data(gene_symbol: str, tf_list: List[str]) -> Tuple[Optional[pd.DataFrame], Optional[str]]:
    """
    Retrieves median tissue expression for target gene and candidate TFs from GTEx Portal API.
    Returns: (DataFrame with columns ['tissue', 'TPM', 'Gene'], error_message_or_None)
    """
    symbols = [gene_symbol.strip()] + [tf.strip() for tf in tf_list if tf.strip()]
    symbols = list(dict.fromkeys(symbols))  # Deduplicate preserving order
    
    gencode_map: Dict[str, str] = {}
    for sym in symbols:
        gid = resolve_gtex_gencode_id(sym)
        if gid:
            gencode_map[sym] = gid
            
    if not gencode_map.get(gene_symbol.strip()):
        return None, f"Could not resolve GENCODE ID for target gene '{gene_symbol}' on GTEx Portal."
        
    expression_dfs = []
    for sym, gid in gencode_map.items():
        expr_data = get_gtex_median_expression(gid)
        if expr_data:
            rows = []
            for entry in expr_data:
                rows.append({
                    "tissue": entry.get("tissueSiteDetailId", "Unknown"),
                    "TPM": float(entry.get("median", 0.0)),
                    "Gene": sym
                })
            if rows:
                expression_dfs.append(pd.DataFrame(rows))
                
    if not expression_dfs:
        return None, "Failed to retrieve median expression data from GTEx API (empty response or offline)."
        
    df_combined = pd.concat(expression_dfs, ignore_index=True)
    return df_combined, None


def compute_tissue_coexpression_summary(
    df_gtex: Optional[pd.DataFrame], 
    gene_symbol: str, 
    tf_list: List[str]
) -> pd.DataFrame:
    """
    Computes cross-tissue Spearman rank correlation of median TPM profiles between
    the target gene and candidate transcription factors across shared GTEx tissue sites.
    
    IMPORTANT SCIENTIFIC DISCLAIMER (DEF-014):
    Cross-tissue median TPM correlation reflects anatomical co-distribution across organ sites,
    NOT individual patient sample co-expression, and does not demonstrate direct causal transcriptional regulation.
    """
    import numpy as np
    import scipy.stats
    
    columns = ["TF", "Target", "Shared_Tissues", "Spearman_Rho", "Spearman_Pval", "Interpretation"]
    if df_gtex is None or df_gtex.empty or "tissue" not in df_gtex.columns or "TPM" not in df_gtex.columns or "Gene" not in df_gtex.columns:
        return pd.DataFrame(columns=columns)
        
    target_clean = gene_symbol.strip()
    target_df = df_gtex[df_gtex["Gene"] == target_clean]
    if target_df.empty:
        return pd.DataFrame(columns=columns)
        
    target_series = target_df.set_index("tissue")["TPM"]
    
    rows = []
    for tf in tf_list:
        tf_clean = tf.strip()
        if not tf_clean or tf_clean == target_clean:
            continue
        tf_df = df_gtex[df_gtex["Gene"] == tf_clean]
        if tf_df.empty:
            continue
        tf_series = tf_df.set_index("tissue")["TPM"]
        
        # Align on shared tissues
        aligned = pd.concat([target_series.rename("Target"), tf_series.rename("TF")], axis=1, join="inner").dropna()
        n_shared = len(aligned)
        
        if n_shared >= 5:
            try:
                rho, pval = scipy.stats.spearmanr(aligned["Target"], aligned["TF"])
                if np.isnan(rho):
                    rho, pval = 0.0, 1.0
                if rho > 0.4 and pval < 0.05:
                    interp = "Tissue Co-Enriched (Anatomical)"
                elif rho < -0.4 and pval < 0.05:
                    interp = "Tissue Mutually Exclusive"
                else:
                    interp = "Tissue Independent"
            except Exception:
                rho, pval, interp = 0.0, 1.0, "Calculation Error"
        else:
            rho, pval, interp = np.nan, np.nan, "Insufficient Shared Tissues (<5)"
            
        rows.append({
            "TF": tf_clean,
            "Target": target_clean,
            "Shared_Tissues": n_shared,
            "Spearman_Rho": round(float(rho), 3) if not np.isnan(rho) else np.nan,
            "Spearman_Pval": round(float(pval), 4) if not np.isnan(pval) else np.nan,
            "Interpretation": interp
        })
        
    return pd.DataFrame(rows, columns=columns)


# ============================================================================
# STRING Database API (Protein-Protein Interactions)
# ============================================================================

def fetch_string_data(
    gene_symbol: str, 
    tf_list: List[str], 
    species: int = 9606
) -> Tuple[Optional[bytes], Optional[pd.DataFrame], Optional[str]]:
    """
    Queries STRING DB REST API v12 for:
    1. Protein ID mapping (resolves symbols to STRING identifiers)
    2. High-resolution protein network visualization (PNG bytes)
    3. Interaction partners with confidence scores (DataFrame)
    
    Returns: (png_image_bytes, partners_dataframe, error_message_or_None)
    """
    identifiers = [gene_symbol.strip()] + [tf.strip() for tf in tf_list if tf.strip()]
    identifiers = list(dict.fromkeys(identifiers))
    
    # 1. Map Identifiers
    map_url = f"{STRING_BASE_URL}/tsv/get_string_ids"
    map_payload = {
        "identifiers": "\r".join(identifiers),
        "species": str(species),
        "caller_identity": "tf-explorer"
    }
    
    mapped_ids = list(identifiers)
    try:
        resp = requests.post(map_url, data=map_payload, timeout=DEFAULT_TIMEOUT)
        if resp.status_code == 200 and resp.text.strip():
            df_map = pd.read_csv(io.StringIO(resp.text), sep="\t")
            if not df_map.empty and "stringId" in df_map.columns:
                mapped = df_map["stringId"].dropna().tolist()
                if mapped:
                    mapped_ids = mapped
    except Exception as e:
        logger.warning(f"STRING ID mapping warning: {e}")
        
    # 2. Network Image
    img_bytes: Optional[bytes] = None
    img_url = f"{STRING_BASE_URL}/highres_image/network"
    img_payload = {
        "identifiers": "\r".join(mapped_ids),
        "species": str(species),
        "caller_identity": "tf-explorer"
    }
    try:
        resp_img = requests.post(img_url, data=img_payload, timeout=15)
        if resp_img.status_code == 200 and resp_img.content:
            img_bytes = resp_img.content
    except Exception as e:
        logger.warning(f"STRING network image fetch failed: {e}")
        
    # 3. Interaction Partners
    df_partners = pd.DataFrame()
    part_url = f"{STRING_BASE_URL}/tsv/interaction_partners"
    part_payload = {
        "identifiers": "\r".join(mapped_ids),
        "species": str(species),
        "limit": "15",
        "caller_identity": "tf-explorer"
    }
    try:
        resp_part = requests.post(part_url, data=part_payload, timeout=DEFAULT_TIMEOUT)
        if resp_part.status_code == 200 and resp_part.text.strip():
            df_partners = pd.read_csv(io.StringIO(resp_part.text), sep="\t")
    except Exception as e:
        logger.warning(f"STRING interaction partners fetch failed: {e}")
        
    if img_bytes is None and df_partners.empty:
        return None, None, "Failed to retrieve STRING network or partner data (network unavailable or no matching interactions)."
        
    return img_bytes, df_partners, None
