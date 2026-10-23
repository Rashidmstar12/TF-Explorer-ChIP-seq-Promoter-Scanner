import requests
import pandas as pd
import logging
from typing import List, Dict, Optional, Tuple, Any
import scipy.stats

# Configure logging
logger = logging.getLogger(__name__)

BASE_URL = "https://www.cbioportal.org/api"

class TCGAClient:
    """
    Client for interacting with cBioPortal API to fetch TCGA expression data.
    """

    # Cache for study metadata to avoid repeated API calls
    _metadata_cache = {}

    @staticmethod
    def fetch_all_tcga_studies():
        """
        Fetches all available TCGA studies from cBioPortal.
        Returns a list of dicts: [{'name': 'Breast Cancer (BRCA)', 'id': 'brca_tcga'}, ...]
        """
        try:
            resp = requests.get("{}/studies".format(BASE_URL), timeout=10)
            if resp.status_code == 200:
                all_studies = resp.json()
                # Filter for TCGA studies (often ending in _tcga or having TCGA in name)
                # Standard convention is id ends with _tcga for PanCancer Atlas or similar
                tcga_studies = []
                for s in all_studies:
                    if s['studyId'].endswith('_tcga'):
                        tcga_studies.append({'name': s['name'], 'id': s['studyId']})
                
                # Sort by name
                tcga_studies.sort(key=lambda x: x['name'])
                return tcga_studies
            return []
        except Exception as e:
            logger.error("Failed to fetch studies: {}".format(e))
            return []

    @staticmethod
    def discover_study_metadata(study_id):
        """
        Dynamically finds the best mRNA molecular profile and sample list for a study.
        Returns (profile_id, sample_list_id) or (None, None).
        """
        if study_id in TCGAClient._metadata_cache:
            return TCGAClient._metadata_cache[study_id]
            
        try:
            # 1. profiles
            p_url = "{}/studies/{}/molecular-profiles".format(BASE_URL, study_id)
            p_resp = requests.get(p_url, timeout=10)
            profile_id = None
            if p_resp.status_code == 200:
                profiles = p_resp.json()
                # Prioritize rna_seq_v2_mrna
                target_p = next((p for p in profiles if p['molecularProfileId'] == "{}_rna_seq_v2_mrna".format(study_id)), None)
                if not target_p:
                    target_p = next((p for p in profiles if "_mrna" in p['molecularProfileId'].lower() and "zscore" not in p['molecularProfileId'].lower()), None)
                
                if target_p:
                    profile_id = target_p['molecularProfileId']
            
            # 2. Sample Lists
            sl_url = "{}/studies/{}/sample-lists".format(BASE_URL, study_id)
            sl_resp = requests.get(sl_url, timeout=10)
            sample_list_id = "{}_all".format(study_id) # Default fallback
            
            if sl_resp.status_code == 200:
                s_lists = sl_resp.json()
                # Try exact match for _all
                target_sl = next((s for s in s_lists if s['sampleListId'] == "{}_all".format(study_id)), None)
                if not target_sl and profile_id:
                     # Try same as profile (some studies use profile ID as list ID)
                    target_sl = next((s for s in s_lists if s['sampleListId'] == profile_id), None)
                if not target_sl:
                    target_sl = next((s for s in s_lists if "_mrna" in s['sampleListId'].lower()), None)
                
                if target_sl:
                    sample_list_id = target_sl['sampleListId']
            
            if profile_id:
                TCGAClient._metadata_cache[study_id] = (profile_id, sample_list_id)
                return profile_id, sample_list_id
            
            return None, None

        except Exception as e:
            logger.error("Metadata discovery failed for {}: {}".format(study_id, e))
            return None, None

    @staticmethod
    def resolve_entrez_ids(genes):
        """
        Resolves gene symbols to Entrez Gene IDs using MyGene.info API.
        Returns a dictionary: {'BRCA1': 672, ...}
        """
        if not genes:
            return {}
            
        try:
            # MyGene.info API
            url = "https://mygene.info/v3/query"
            headers = {'content-type': 'application/x-www-form-urlencoded'}
            params = {
                'q': ",".join(genes),
                'scopes': 'symbol,alias', # Check symbol and aliases
                'fields': 'entrezgene,symbol',
                'species': 'human'
            }
            
            resp = requests.post(url, data=params, headers=headers, timeout=10)
            
            if resp.status_code == 200:
                results = resp.json()
                gene_map = {}
                for res in results:
                    # Check if 'entrezgene' field exists (it might be missing if not found)
                    if 'entrezgene' in res:
                        # Use the query term as the key to map back to input
                        query_term = res.get('query')
                        # Or use the returned symbol if you prefer canonical
                        canonical_symbol = res.get('symbol', query_term)
                        
                        # Map input symbol to ID. 
                        # Note: One input might match multiple or none. 
                        # This simple logic maps query -> ID
                        if query_term:
                           gene_map[query_term] = int(res['entrezgene'])
                           
                        # Also map canonical symbol just in case input was an alias
                        if canonical_symbol and canonical_symbol != query_term:
                            gene_map[canonical_symbol] = int(res['entrezgene'])
                            
                return gene_map
            else:
                logger.error("MyGene info query failed: {}".format(resp.status_code))
                return {}
                
        except Exception as e:
            logger.error("Failed to resolve gene IDs: {}".format(e))
            return {}

    @staticmethod
    def fetch_expression_with_status(genes: List[str], study_id: str) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """
        Fetches mRNA expression data for the given genes in the specified study with explicit status reporting.
        Never falls back to synthetic or fabricated patient distributions (DEF-023).
        Returns: (df_pivot, status_dict)
        """
        if not genes or not study_id:
            return pd.DataFrame(), {
                "status": "UNAVAILABLE",
                "message": "Gene list or TCGA study ID was empty.",
                "study_id": study_id,
                "genes_queried": genes or []
            }

        # Discover Metadata
        profile_id, sample_list_id = TCGAClient.discover_study_metadata(study_id)
        if not profile_id:
            msg = f"Could not find valid mRNA expression profile for TCGA study '{study_id}'."
            logger.warning(msg)
            return pd.DataFrame(), {
                "status": "UNAVAILABLE",
                "message": msg,
                "study_id": study_id,
                "genes_queried": genes
            }
        
        # 1. Resolve Genes
        gene_map = TCGAClient.resolve_entrez_ids(genes)
        if not gene_map:
            msg = f"No Entrez IDs resolved for queried genes: {genes}."
            logger.warning(msg)
            return pd.DataFrame(), {
                "status": "UNAVAILABLE",
                "message": msg,
                "study_id": study_id,
                "genes_queried": genes
            }
        
        entrez_ids = list(gene_map.values())
        
        # 2. Fetch Data
        url = f"{BASE_URL}/molecular-profiles/{profile_id}/molecular-data/fetch"
        payload = {
            "entrezGeneIds": entrez_ids,
            "sampleListId": sample_list_id
        }
        
        logger.info(f"Fetching TCGA data for {genes} from {profile_id} (list: {sample_list_id})...")
        
        try:
            resp = requests.post(url, json=payload, timeout=30)
            if resp.status_code != 200:
                msg = f"cBioPortal API request failed with HTTP status {resp.status_code}."
                logger.error(msg)
                return pd.DataFrame(), {
                    "status": "UNAVAILABLE",
                    "message": msg,
                    "study_id": study_id,
                    "genes_queried": genes
                }
            data = resp.json()
            
            # 3. Parse Response
            records = []
            id_to_symbol = {v: k for k, v in gene_map.items()}
            
            for item in data:
                gid = item.get('entrezGeneId')
                symbol = id_to_symbol.get(gid, f"Entrez_{gid}")
                records.append({
                    "sampleId": item.get('sampleId'),
                    "gene": symbol,
                    "expression": item.get('value')
                })
                
            if not records:
                msg = f"No expression records returned from cBioPortal for {genes} in '{study_id}'."
                return pd.DataFrame(), {
                    "status": "UNAVAILABLE",
                    "message": msg,
                    "study_id": study_id,
                    "genes_queried": genes
                }
                
            df_long = pd.DataFrame(records)
            df_pivot = df_long.pivot(index='sampleId', columns='gene', values='expression').dropna()
            
            return df_pivot, {
                "status": "SUCCESS",
                "message": f"Successfully retrieved {len(df_pivot)} patient sample profiles.",
                "study_id": study_id,
                "genes_queried": genes,
                "samples_count": len(df_pivot)
            }
            
        except requests.exceptions.RequestException as e:
            msg = f"Network error fetching TCGA expression data: {e}"
            logger.error(msg)
            return pd.DataFrame(), {
                "status": "UNAVAILABLE",
                "message": msg,
                "study_id": study_id,
                "genes_queried": genes
            }

    @staticmethod
    def fetch_expression_data(genes, study_id):
        """
        Fetches mRNA expression data for the given genes in the specified study.
        Returns pd.DataFrame (empty if unavailable).
        """
        df, _ = TCGAClient.fetch_expression_with_status(genes, study_id)
        return df

    @staticmethod
    def calculate_correlation(df: pd.DataFrame, tf: str, target: str) -> Dict[str, Any]:
        """
        Calculates Pearson and Spearman correlations between TF and Target.
        Returns explicit status dictionary with no synthetic fallback (DEF-023).
        """
        if df is None or df.empty or tf not in df.columns or target not in df.columns:
            return {
                "status": "UNAVAILABLE",
                "message": f"Cannot compute correlation: required columns ('{tf}', '{target}') missing or dataset empty.",
                "pearson_r": None,
                "pearson_p": None,
                "spearman_r": None,
                "spearman_p": None,
                "n_samples": 0
            }
            
        x = df[tf]
        y = df[target]
        
        if len(x) < 3:
            return {
                "status": "UNAVAILABLE",
                "message": f"Insufficient sample size (N={len(x)}) to compute meaningful correlation.",
                "pearson_r": None,
                "pearson_p": None,
                "spearman_r": None,
                "spearman_p": None,
                "n_samples": len(x)
            }
        
        try:
            pearson_r, pearson_p = scipy.stats.pearsonr(x, y)
            spearman_r, spearman_p = scipy.stats.spearmanr(x, y)
            
            return {
                "status": "SUCCESS",
                "message": "Correlation computed successfully.",
                "pearson_r": float(pearson_r),
                "pearson_p": float(pearson_p),
                "spearman_r": float(spearman_r),
                "spearman_p": float(spearman_p),
                "n_samples": int(len(df))
            }
        except Exception as e:
            msg = f"Correlation calculation error: {e}"
            logger.error(msg)
            return {
                "status": "UNAVAILABLE",
                "message": msg,
                "pearson_r": None,
                "pearson_p": None,
                "spearman_r": None,
                "spearman_p": None,
                "n_samples": len(df)
            }
