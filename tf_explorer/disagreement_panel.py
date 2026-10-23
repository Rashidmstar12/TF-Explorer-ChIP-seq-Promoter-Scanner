"""
Streamlit Disagreement Explorer Panel for TF-Explorer Phase 4B-R.

Provides an interactive workbench for exploring, diagnosing, and exporting
evidence disagreement across transcripts, promoter windows, experiments, and contexts.
All presets and custom analyses pass strictly through build_validated_configuration
for fail-closed primary BED authentication and verifiable peak overlap calculations.
"""

import json
from pathlib import Path
from typing import Dict, Any, List, Optional
import streamlit as st
import pandas as pd

from tf_explorer.disagreement_taxonomy import (
    EvidenceConfiguration,
    EvidenceInterpretation,
    DisagreementClassification
)
from tf_explorer.disagreement_factory import build_validated_configuration
from tf_explorer.disagreement_matrix import DisagreementMatrix
from tf_explorer.peak_authenticator import PrimaryPeakAuthenticator
from tf_explorer.disagreement_provenance import AUTHORITATIVE_ENCODE_CATALOG


def get_file_metadata_from_catalog(file_acc: str) -> Dict[str, str]:
    """Generates researcher-facing metadata directly from the authoritative verified catalog."""
    entry = AUTHORITATIVE_ENCODE_CATALOG.get(file_acc)
    if not entry:
        raise KeyError(f"File accession {file_acc} not found in authoritative ENCODE catalog")
    return {
        "tf": entry["target_tf"],
        "exp": entry["experiment_accession"],
        "file": entry["file_accession"],
        "bio": entry["biosample"],
        "treat": entry["treatment_context"],
        "mod": entry.get("genetic_modification_context", "None")
    }


PRESET_CONFIGURATIONS = {
    "CDK4 / E2F1 (MCF-7 vs HeLa-S3 Biological Comparison)": {
        "gene": "CDK4",
        "tf": "E2F1",
        "desc": "Reconciled comparison between MCF-7 (HA-tagged) and HeLa-S3 (HA-tagged). Both exhibit genuine core promoter binding (~302-314 bp).",
        "configs": [
            {
                "config_id": "CDK4_MCF7_Core",
                "transcript_id": "ENST00000257904.11",
                "tss": 57752310, "strand": "-", "is_mane": True,
                "window": "Core [-500, +100]", "w_start": -500, "w_end": 100,
                "exp": "ENCSR000EWX", "file": "ENCFF115KRO",
                "bio": "MCF-7", "treat": "Untreated", "mod": "HA-tagged", "chrom": "chr12"
            },
            {
                "config_id": "CDK4_HeLa_Core",
                "transcript_id": "ENST00000257904.11",
                "tss": 57752310, "strand": "-", "is_mane": True,
                "window": "Core [-500, +100]", "w_start": -500, "w_end": 100,
                "exp": "ENCSR000EVM", "file": "ENCFF004ZKN",
                "bio": "HeLa-S3", "treat": "Untreated", "mod": "HA-tagged", "chrom": "chr12"
            }
        ]
    },
    "TFF1 / ESR1 (Promoter Window Arithmetic & Flank Sensitivity)": {
        "gene": "TFF1",
        "tf": "ESR1",
        "desc": "Demonstrates window flank arithmetic: standard [-2000, +500] covers full 691 bp peak; downstream-only [0, +500] retains 140 bp partial overlap.",
        "configs": [
            {
                "config_id": "TFF1_Standard_Window",
                "transcript_id": "ENST00000291527.3",
                "tss": 42366535, "strand": "-", "is_mane": True,
                "window": "Standard [-2000, +500]", "w_start": -2000, "w_end": 500,
                "exp": "ENCSR463GOT", "file": "ENCFF185UUT",
                "bio": "MCF-7", "treat": "Unknown / Unspecified", "mod": "eGFP-tagged (C-terminal)", "chrom": "chr21"
            },
            {
                "config_id": "TFF1_Downstream_Only",
                "transcript_id": "ENST00000291527.3",
                "tss": 42366535, "strand": "-", "is_mane": True,
                "window": "Downstream Only [0, +500]", "w_start": 0, "w_end": 500,
                "exp": "ENCSR463GOT", "file": "ENCFF185UUT",
                "bio": "MCF-7", "treat": "Unknown / Unspecified", "mod": "eGFP-tagged (C-terminal)", "chrom": "chr21"
            }
        ]
    },
    "TP53 / E2F1 (Alternative Transcript TSS Sensitivity)": {
        "gene": "TP53",
        "tf": "E2F1",
        "desc": "Canonical MANE Select TSS (7,687,550) has 336 bp peak overlap, whereas alternative TSS (7,669,402) has 0 bp overlap.",
        "configs": [
            {
                "config_id": "TP53_MANE_Select",
                "transcript_id": "ENST00000269305.9",
                "tss": 7687550, "strand": "-", "is_mane": True,
                "window": "Core [-500, +100]", "w_start": -500, "w_end": 100,
                "exp": "ENCSR000EWX", "file": "ENCFF115KRO",
                "bio": "MCF-7", "treat": "Untreated", "mod": "HA-tagged", "chrom": "chr17"
            },
            {
                "config_id": "TP53_Alternative_Transcript",
                "transcript_id": "ENST00000455683.1",
                "tss": 7669402, "strand": "-", "is_mane": False,
                "window": "Core [-500, +100]", "w_start": -500, "w_end": 100,
                "exp": "ENCSR000EWX", "file": "ENCFF115KRO",
                "bio": "MCF-7", "treat": "Untreated", "mod": "HA-tagged", "chrom": "chr17"
            }
        ]
    },
    "TP53 / CTCF (Promoter Window Boundary Sensitivity)": {
        "gene": "TP53",
        "tf": "CTCF",
        "desc": "CTCF peak at -644 bp from TSS is included in extended window [-2000, +500] (54 bp overlap) but excluded in core window [-500, +100] (0 bp overlap).",
        "configs": [
            {
                "config_id": "TP53_CTCF_Extended",
                "transcript_id": "ENST00000269305.9",
                "tss": 7687550, "strand": "-", "is_mane": True,
                "window": "Extended [-2000, +500]", "w_start": -2000, "w_end": 500,
                "exp": "ENCSR000DMS", "file": "ENCFF685HMV",
                "bio": "MCF-7", "treat": "17β-estradiol (100 nM, 45 min)", "mod": "None", "chrom": "chr17"
            },
            {
                "config_id": "TP53_CTCF_Core",
                "transcript_id": "ENST00000269305.9",
                "tss": 7687550, "strand": "-", "is_mane": True,
                "window": "Core [-500, +100]", "w_start": -500, "w_end": 100,
                "exp": "ENCSR000DMS", "file": "ENCFF685HMV",
                "bio": "MCF-7", "treat": "17β-estradiol (100 nM, 45 min)", "mod": "None", "chrom": "chr17"
            }
        ]
    },
    "CDKN1A / E2F1 (Unseen Locus Multi-Window Flank Analysis)": {
        "gene": "CDKN1A",
        "tf": "E2F1",
        "desc": "Accurate flank arithmetic on positive-strand CDKN1A: Proximal [-1000, +200] (425 bp) -> Downstream [+50, +500] (54 bp) -> Distal [+600, +1500] (0 bp).",
        "configs": [
            {
                "config_id": "CDKN1A_Proximal",
                "transcript_id": "ENST00000244741.10",
                "tss": 36678714, "strand": "+", "is_mane": True,
                "window": "Proximal [-1000, +200]", "w_start": -1000, "w_end": 200,
                "exp": "ENCSR000EWX", "file": "ENCFF115KRO",
                "bio": "MCF-7", "treat": "Untreated", "mod": "HA-tagged", "chrom": "chr6"
            },
            {
                "config_id": "CDKN1A_Downstream_50",
                "transcript_id": "ENST00000244741.10",
                "tss": 36678714, "strand": "+", "is_mane": True,
                "window": "Downstream [+50, +500]", "w_start": 50, "w_end": 500,
                "exp": "ENCSR000EWX", "file": "ENCFF115KRO",
                "bio": "MCF-7", "treat": "Untreated", "mod": "HA-tagged", "chrom": "chr6"
            },
            {
                "config_id": "CDKN1A_Distal",
                "transcript_id": "ENST00000244741.10",
                "tss": 36678714, "strand": "+", "is_mane": True,
                "window": "Distal [+600, +1500]", "w_start": 600, "w_end": 1500,
                "exp": "ENCSR000EWX", "file": "ENCFF115KRO",
                "bio": "MCF-7", "treat": "Untreated", "mod": "HA-tagged", "chrom": "chr6"
            }
        ]
    },
    "Custom Investigation": None
}

SUPPORTED_GENES = {
    "CDK4": {
        "chrom": "chr12", "strand": "-", "canonical_tss": 57752310,
        "transcripts": ["ENST00000257904.11 (MANE Select)", "ENST00000544414.1"],
        "available_files": [
            get_file_metadata_from_catalog("ENCFF115KRO"),
            get_file_metadata_from_catalog("ENCFF004ZKN"),
            get_file_metadata_from_catalog("ENCFF142BHB"),
            get_file_metadata_from_catalog("ENCFF030DLI"),
        ]
    },
    "TFF1": {
        "chrom": "chr21", "strand": "-", "canonical_tss": 42366535,
        "transcripts": ["ENST00000291527.3 (MANE Select)"],
        "available_files": [
            get_file_metadata_from_catalog("ENCFF185UUT"),
            get_file_metadata_from_catalog("ENCFF142BHB"),
            get_file_metadata_from_catalog("ENCFF030DLI"),
        ]
    },
    "TP53": {
        "chrom": "chr17", "strand": "-", "canonical_tss": 7687550,
        "transcripts": ["ENST00000269305.9 (MANE Select)", "ENST00000455683.1 (TSS 7,669,402)"],
        "available_files": [
            get_file_metadata_from_catalog("ENCFF115KRO"),
            get_file_metadata_from_catalog("ENCFF685HMV"),
            get_file_metadata_from_catalog("ENCFF142BHB"),
        ]
    },
    "CDKN1A": {
        "chrom": "chr6", "strand": "+", "canonical_tss": 36678714,
        "transcripts": ["ENST00000244741.10 (MANE Select)"],
        "available_files": [
            get_file_metadata_from_catalog("ENCFF115KRO"),
            get_file_metadata_from_catalog("ENCFF142BHB"),
        ]
    }
}


def run_disagreement_explorer_panel():
    st.subheader("⚖️ Regulatory Evidence Disagreement Explorer")
    st.markdown(
        "Diagnose and explain conflicting promoter-associated transcription factor evidence across "
        "transcripts, promoter windows, experiments, biological contexts, and processing settings. "
        "**Strictly fail-closed: All configurations derive evidence states from authenticated primary BED bytes.**"
    )

    st.markdown("### Step 1: Select Preset Benchmark Case or Custom Investigation")
    preset_choice = st.selectbox("Benchmark Presets", list(PRESET_CONFIGURATIONS.keys()))
    preset_data = PRESET_CONFIGURATIONS[preset_choice]

    authenticator = PrimaryPeakAuthenticator()

    if preset_data is not None:
        st.info(f"**Case Description**: {preset_data['desc']}")
        gene_val = preset_data["gene"]
        tf_val = preset_data["tf"]
        matrix = DisagreementMatrix(gene_symbol=gene_val, tf_target=tf_val)

        for raw_cfg in preset_data["configs"]:
            cfg = build_validated_configuration(
                gene_symbol=gene_val,
                transcript_id=raw_cfg["transcript_id"],
                canonical_tss=raw_cfg["tss"],
                strand=raw_cfg["strand"],
                tf_target=tf_val,
                encode_experiment_accession=raw_cfg["exp"],
                encode_file_accession=raw_cfg["file"],
                biosample=raw_cfg["bio"],
                treatment_context=raw_cfg["treat"],
                genetic_modification_context=raw_cfg["mod"],
                window_name=raw_cfg["window"],
                window_start_offset=raw_cfg["w_start"],
                window_end_offset=raw_cfg["w_end"],
                chromosome=raw_cfg.get("chrom", ""),
                config_id=raw_cfg["config_id"],
                authenticator=authenticator,
                is_mane_select=raw_cfg["is_mane"],
                mane_status="MANE_SELECT_VERIFIED" if raw_cfg["is_mane"] else "MANE_UNVERIFIED"
            )
            matrix.add_configuration(cfg)

    else:
        # Custom Investigation Mode
        st.markdown("#### 🔬 Custom Evidence Configuration Workbench")
        st.caption("Demonstration Catalog (4 benchmark genes: CDK4, TFF1, TP53, CDKN1A). Strict fail-closed verification against authoritative ENCODE primary BED files.")
        st.markdown("Construct and evaluate an authenticated evidence matrix for supported benchmark loci.")

        col_g1, col_g2 = st.columns(2)
        with col_g1:
            gene_sel = st.selectbox("Target Gene", list(SUPPORTED_GENES.keys()))
            gene_info = SUPPORTED_GENES[gene_sel]
        with col_g2:
            tx_sel_str = st.selectbox("Transcript", gene_info["transcripts"])
            tx_id = tx_sel_str.split()[0]
            is_mane = "MANE Select" in tx_sel_str
            tss_val = 7669402 if "7,669,402" in tx_sel_str else gene_info["canonical_tss"]

        tf_options = sorted(list({f["tf"] for f in gene_info["available_files"]}))
        col_t1, col_t2 = st.columns(2)
        with col_t1:
            tf_sel = st.selectbox("Transcription Factor Target", tf_options)
        with col_t2:
            matched_files = [f for f in gene_info["available_files"] if f["tf"] == tf_sel]
            file_options = [f"{f['file']} ({f['exp']} | {f['bio']} | {f['treat']} | {f.get('mod', 'None')})" for f in matched_files]
            file_sels = st.multiselect("Select ENCODE Files to Compare", file_options, default=file_options[:2] if len(file_options) >= 2 else file_options)

        col_w1, col_w2 = st.columns(2)
        with col_w1:
            w_start = st.number_input("Window 5' Offset (relative to TSS)", value=-1000, step=100)
        with col_w2:
            w_end = st.number_input("Window 3' Offset (relative to TSS)", value=200, step=50)

        matrix = DisagreementMatrix(gene_symbol=gene_sel, tf_target=tf_sel)

        for sel_str in file_sels:
            acc = sel_str.split()[0]
            f_meta = next(f for f in matched_files if f["file"] == acc)
            cfg = build_validated_configuration(
                gene_symbol=gene_sel,
                transcript_id=tx_id,
                canonical_tss=tss_val,
                strand=gene_info["strand"],
                tf_target=tf_sel,
                encode_experiment_accession=f_meta["exp"],
                encode_file_accession=acc,
                biosample=f_meta["bio"],
                treatment_context=f_meta["treat"],
                genetic_modification_context=f_meta.get("mod", "None"),
                window_name=f"Window [{w_start}, {w_end}]",
                window_start_offset=int(w_start),
                window_end_offset=int(w_end),
                chromosome=gene_info["chrom"],
                config_id=f"{gene_sel}_{f_meta['bio']}_{acc}",
                authenticator=authenticator,
                is_mane_select=is_mane,
                mane_status="MANE_SELECT_VERIFIED" if is_mane else "MANE_UNVERIFIED"
            )
            matrix.add_configuration(cfg)

    # Step 2: Matrix Execution and Display
    st.markdown("### Step 2: Evidence Matrix & Cause-of-Disagreement Analysis")
    stats = matrix.get_summary_statistics()
    analysis = matrix.analyze()

    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("Configurations", stats["total_configurations"])
    col2.metric("Independent Exps", stats["unique_experiments_count"])
    col3.metric("Positive Cells", stats["positive_configurations_count"])
    col4.metric("Negative Cells", stats["negative_configurations_count"])
    col5.metric("Excluded Cells", stats["excluded_configurations_count"])

    st.info(f"**Scientific Interpretation**: {stats['scientific_interpretation']}")

    st.markdown("### Configuration Evidence Grid")
    df_grid = matrix.to_dataframe()
    st.dataframe(df_grid, use_container_width=True)

    if analysis["diagnoses"]:
        st.markdown("### Detailed Cause-of-Disagreement Diagnoses")
        for d in analysis["diagnoses"]:
            with st.expander(f"Disagreement: `{d['config_a_id']}` vs `{d['config_b_id']}`", expanded=True):
                st.write(f"**Primary Causes**: {', '.join(d['primary_causes'])}")
                st.write(f"**Smallest Observed Difference**: {d['smallest_observed_difference']}")
                st.write(f"**Explanation**: {d['scientific_explanation']}")
                if d["causal_attribution_caveat"]:
                    st.warning(f"**Attribution Caveat**: {d['causal_attribution_caveat']}")

    st.markdown("### Step 3: Reproducibility Exports")
    col_e1, col_e2, col_e3 = st.columns(3)
    gene_export = matrix.gene_symbol
    tf_export = matrix.tf_target
    with col_e1:
        csv_data = df_grid.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download Matrix CSV",
            csv_data,
            file_name=f"{gene_export}_{tf_export}_disagreement_matrix.csv",
            mime="text/csv"
        )
    with col_e2:
        json_data = json.dumps(analysis, indent=2).encode("utf-8")
        st.download_button(
            "Download Diagnosis JSON",
            json_data,
            file_name=f"{gene_export}_{tf_export}_disagreement_analysis.json",
            mime="application/json"
        )
    with col_e3:
        md_data = matrix.to_markdown().encode("utf-8")
        st.download_button(
            "Download Markdown Dossier",
            md_data,
            file_name=f"{gene_export}_{tf_export}_disagreement_dossier.md",
            mime="text/markdown"
        )
