import requests
import os
import logging
from typing import List, Dict, Optional

# Configure logging
logger = logging.getLogger(__name__)

ENCODE_BASE_URL = "https://www.encodeproject.org"

TF_ALIASES = {
    "CREB": "CREB1",
    "P53": "TP53",
    "C-MYC": "MYC"
}

def search_encode(target_name: str, assay_type: str = "TF ChIP-seq", organism: str = "Homo sapiens") -> List[Dict]:
    """
    General ENCODE search supporting TF ChIP-seq, Histone ChIP-seq, and DNase-seq.
    """
    logger.info(f"Searching ENCODE for {target_name} ({assay_type}, {organism})...")
    
    # 1. Base parameters
    params = {
        "type": "Experiment",
        "assay_title": assay_type,
        "replicates.library.biosample.donor.organism.scientific_name": organism,
        "status": "released",
        "limit": "all",
        "field": ["accession", "files", "biosample_ontology", "description"]
    }
    
    # 2. Add target for targeted assays
    if assay_type in ["TF ChIP-seq", "Histone ChIP-seq"]:
        # Resolve alias if TF ChIP-seq
        if assay_type == "TF ChIP-seq":
            original_name = target_name
            target_name = TF_ALIASES.get(target_name.upper(), target_name)
            if target_name != original_name:
                logger.info(f"Resolved alias: {original_name} -> {target_name}")
        params["target.label"] = target_name
    
    headers = {"Accept": "application/json"}
    
    try:
        response = requests.get(f"{ENCODE_BASE_URL}/search/", params=params, headers=headers)
        response.raise_for_status()
        data = response.json()
        
        experiments = data.get("@graph", [])
        logger.info(f"Found {len(experiments)} experiments for {target_name} ({assay_type}).")
        
        results = []
        for exp in experiments:
            if not isinstance(exp, dict):
                continue
            
            biosample = exp.get("biosample_ontology", {}).get("term_name", "Unknown")
            files = exp.get("files", [])
            for f in files:
                if not isinstance(f, dict):
                    continue
                
                # Setup assay-specific file filters
                if assay_type == "TF ChIP-seq":
                    valid_output_types = [
                        "optimal idr thresholded peaks", 
                        "IDR thresholded peaks", 
                        "conservative IDR thresholded peaks"
                    ]
                    valid_file_format_types = ["narrowPeak"]
                elif assay_type == "Histone ChIP-seq":
                    valid_output_types = [
                        "reproducible peaks",
                        "peaks",
                        "optimal idr thresholded peaks", 
                        "IDR thresholded peaks", 
                        "conservative IDR thresholded peaks"
                    ]
                    valid_file_format_types = ["narrowPeak", "broadPeak"]
                elif assay_type == "DNase-seq":
                    valid_output_types = ["peaks", "hotspots"]
                    valid_file_format_types = ["narrowPeak", "broadPeak"]
                else:
                    valid_output_types = ["peaks"]
                    valid_file_format_types = ["narrowPeak", "broadPeak"]
                
                output_type = f.get("output_type", "")
                file_format = f.get("file_format", "")
                file_format_type = f.get("file_format_type", "")
                assembly = f.get("assembly", "")
                status = f.get("status", "")
                
                if (output_type in valid_output_types and
                    file_format == "bed" and
                    file_format_type in valid_file_format_types and
                    assembly == "GRCh38" and
                    status == "released"):
                    
                    results.append({
                        "file_accession": f.get("accession"),
                        "download_url": f"{ENCODE_BASE_URL}{f.get('href')}",
                        "dataset_accession": exp.get("accession"),
                        "assembly": assembly,
                        "biosample": biosample,
                        "description": exp.get("description", "No description")
                    })
        
        logger.info(f"Found {len(results)} peak files for {target_name} from {len(experiments)} experiments.")
        return results
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Error querying ENCODE: {e}")
        return []

def search_encode_tf_chipseq(tf_name: str, organism: str = "Homo sapiens") -> List[Dict]:
    """Backward-compatible wrapper for TF ChIP-seq."""
    return search_encode(tf_name, "TF ChIP-seq", organism)

def download_peak_file(url: str, local_path: str) -> Optional[str]:
    """
    Downloads a file from a URL to a local path with retries.
    Returns the local path if successful, None otherwise.
    """
    import time
    max_retries = 3
    for attempt in range(max_retries):
        try:
            logger.info(f"Downloading {url} to {local_path} (Attempt {attempt+1}/{max_retries})...")
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()
            
            with open(local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
                    
            return local_path
        except requests.exceptions.RequestException as e:
            logger.warning(f"Failed to download {url}: {e}")
            if attempt < max_retries - 1:
                time.sleep(2 * (attempt + 1)) # Exponential backoff
            else:
                logger.error(f"Given up on {url} after {max_retries} attempts.")
                return None

def search_region(chrom: str, start: int, end: int, assembly: str = "GRCh38") -> List[Dict]:
    """
    Searches ENCODE for experiments/files overlapping a specific genomic region.
    Returns a list of unique TFs found in that region with metadata.
    """
    import time
    
    # Format region string
    # ENCODE expects 'chrN:start-end'
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
        
    region_str = f"{chrom}:{start}-{end}"
    logger.info(f"Searching ENCODE region: {region_str} ({assembly})...")
    
    url = f"{ENCODE_BASE_URL}/search/"
    params = {
        "type": "Experiment",
        "assay_title": "TF ChIP-seq",
        "region": region_str,
        "genome": assembly,
        "format": "json",
        "limit": "all",
        "status": "released"
    }
    
    headers = {
        "Accept": "application/json",
        "User-Agent": "TF-Explorer/2.0 (research-tool)"
    }
    
    max_retries = 3
    for attempt in range(max_retries):
        try:
            response = requests.get(url, params=params, headers=headers, timeout=60)
            response.raise_for_status()
            data = response.json()
            
            items = data.get("@graph", [])
            logger.info(f"Found {len(items)} experiments in region search.")
            
            # Aggregate by Target (TF)
            tf_stats = {}
            
            for exp in items:
                # Items are Experiments now
                target = exp.get('target')
                if not target or not isinstance(target, dict):
                    continue
                    
                tf_name = target.get('label')
                if not tf_name:
                    continue
                    
                biosample = exp.get("biosample_ontology", {}).get("term_name", "Unknown")
                
                if tf_name not in tf_stats:
                    tf_stats[tf_name] = {
                        'tf_name': tf_name,
                        'count': 0,
                        'biosamples': set(),
                        'experiments': set()
                    }
                
                # Count files? We don't have file count directly here without fetching each exp.
                # But we can just count experiments.
                # Or we can look at 'files' list if included in response (it usually is for search).
                # Let's just count experiments for now as it's faster.
                tf_stats[tf_name]['count'] += 1 # Increment experiment count as proxy
                tf_stats[tf_name]['biosamples'].add(biosample)
                tf_stats[tf_name]['experiments'].add(exp['accession'])

            results = []
            for tf, stats in tf_stats.items():
                results.append({
                    'tf_name': tf,
                    'file_count': len(stats['experiments']), # Use exp count
                    'experiment_count': len(stats['experiments']),
                    'biosamples': sorted(list(stats['biosamples']))
                })
                
            results.sort(key=lambda x: x['experiment_count'], reverse=True)
            return results
            
        except requests.exceptions.RequestException as e:
            logger.warning(f"Attempt {attempt+1}/{max_retries} failed for region search: {e}")
            if attempt < max_retries - 1:
                time.sleep(2 * (attempt + 1))
            else:
                logger.error(f"Error querying ENCODE region search after retries: {e}")
                return None
