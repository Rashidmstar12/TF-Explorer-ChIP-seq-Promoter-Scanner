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

def search_encode_tf_chipseq(tf_name: str, organism: str = "Homo sapiens") -> List[Dict]:
    """
    Searches ENCODE for TF ChIP-seq experiments for a given TF and organism.
    Filters for 'optimal idr thresholded peaks' and 'bed narrowPeak' file type.
    """
    # Resolve alias
    original_name = tf_name
    tf_name = TF_ALIASES.get(tf_name.upper(), tf_name)
    
    if tf_name != original_name:
        logger.info(f"Resolved alias: {original_name} -> {tf_name}")

    logger.info(f"Searching ENCODE for {tf_name} ({organism})...")
    
    # 1. Search for Experiments
    params = {
        "type": "Experiment",
        "assay_title": "TF ChIP-seq",
        "target.label": tf_name,
        "replicates.library.biosample.donor.organism.scientific_name": organism,
        "status": "released",
        "limit": "all",
        "field": ["accession", "files", "biosample_ontology", "description"]
    }
    
    headers = {"Accept": "application/json"}
    
    try:
        response = requests.get(f"{ENCODE_BASE_URL}/search/", params=params, headers=headers)
        response.raise_for_status()
        data = response.json()
        
        experiments = data.get("@graph", [])
        logger.info(f"Found {len(experiments)} experiments for {tf_name}.")
        
        results = []
        for exp in experiments:
            if not isinstance(exp, dict):
                logger.warning(f"Experiment entry is not a dict: {exp}")
                continue
            
            # Get biosample from experiment
            biosample = exp.get("biosample_ontology", {}).get("term_name", "Unknown")
                
            files = exp.get("files", [])
            for f in files:
                if not isinstance(f, dict):
                    continue
                    
                # Filter for optimal/IDR thresholded peaks, bed, narrowPeak, GRCh38
                # Relaxed to accept "IDR thresholded peaks" as well, as "optimal" might be archived.
                valid_output_types = [
                    "optimal idr thresholded peaks", 
                    "IDR thresholded peaks", 
                    "conservative IDR thresholded peaks"
                ]
                
                if (f.get("output_type") in valid_output_types and
                    f.get("file_format") == "bed" and
                    f.get("file_format_type") == "narrowPeak" and
                    f.get("assembly") == "GRCh38" and
                    f.get("status") == "released"):
                    
                    results.append({
                        "file_accession": f.get("accession"),
                        "download_url": f"{ENCODE_BASE_URL}{f.get('href')}",
                        "dataset_accession": exp.get("accession"),
                        "assembly": f.get("assembly"),
                        "biosample": biosample,
                        "description": exp.get("description", "No description")
                    })
        
        logger.info(f"Found {len(results)} peak files for {tf_name} from {len(experiments)} experiments.")
        return results
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Error querying ENCODE: {e}")
        return []

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
