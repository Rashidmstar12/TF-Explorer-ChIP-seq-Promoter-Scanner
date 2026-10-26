"""
Multidimensional Regulatory Evidence Matrix for TF-Explorer Phase 4B-R.

Constructs, manages, and exports the configuration-level matrix across:
Transcript x Promoter Window x Experiment x Biological Context.
"""

from typing import List, Dict, Any, Optional, Set, Tuple
from dataclasses import dataclass, field
import json
import pandas as pd

from tf_explorer.disagreement_taxonomy import (
    EvidenceConfiguration,
    EvidenceInterpretation,
    DisagreementClassification
)
from tf_explorer.disagreement_engine import DisagreementEngine


class DisagreementMatrix:
    """
    Multidimensional evidence matrix spanning analytical choices and experimental conditions.
    Decouples configuration space from independent biological experiments to prevent multiverse fallacy.
    """

    def __init__(self, gene_symbol: str, tf_target: str, allow_synthetic_fixtures: bool = False):
        self.gene_symbol = gene_symbol.upper()
        self.tf_target = tf_target.upper()
        self.configurations: List[EvidenceConfiguration] = []
        self.allow_synthetic_fixtures = allow_synthetic_fixtures
        self.engine = DisagreementEngine(allow_synthetic_fixtures=allow_synthetic_fixtures)
        self._analysis_cache: Optional[Dict[str, Any]] = None

    def add_configuration(self, config: EvidenceConfiguration) -> None:
        """Adds an evaluated configuration to the matrix."""
        self.configurations.append(config)
        self._analysis_cache = None

    def clear(self) -> None:
        """Clears all configurations."""
        self.configurations.clear()
        self._analysis_cache = None

    def analyze(self) -> Dict[str, Any]:
        """Runs the cause-of-disagreement engine across all configurations in the matrix."""
        if self._analysis_cache is None:
            self._analysis_cache = self.engine.evaluate_set(self.configurations, allow_synthetic_fixtures=self.allow_synthetic_fixtures)
        return self._analysis_cache

    def get_summary_statistics(self) -> Dict[str, Any]:
        """
        Returns structured summary metrics.
        Explicitly distinguishes configuration counts from independent experiment counts.
        """
        analysis = self.analyze()
        unique_exps = {c.encode_experiment_accession for c in self.configurations if c.encode_experiment_accession}
        unique_txs = {c.transcript_id for c in self.configurations if c.transcript_id}
        unique_windows = {c.window_name for c in self.configurations if c.window_name}
        unique_biosamples = {c.biosample for c in self.configurations if c.biosample}

        return {
            "gene_symbol": self.gene_symbol,
            "tf_target": self.tf_target,
            "total_configurations": len(self.configurations),
            "unique_transcripts_count": len(unique_txs),
            "unique_windows_count": len(unique_windows),
            "unique_experiments_count": len(unique_exps),
            "unique_biosamples_count": len(unique_biosamples),
            "positive_configurations_count": analysis["positive_configurations_count"],
            "negative_configurations_count": analysis["negative_configurations_count"],
            "excluded_configurations_count": analysis["excluded_configurations_count"],
            "overall_status": analysis["overall_status"],
            "observed_causes": analysis["observed_causes"],
            "scientific_interpretation": analysis["scientific_interpretation"]
        }

    def to_dataframe(self) -> pd.DataFrame:
        """Exports configurations to a Pandas DataFrame."""
        rows = []
        for c in self.configurations:
            status_val = c.evidence_status.value if isinstance(c.evidence_status, EvidenceInterpretation) else str(c.evidence_status)
            rows.append({
                "Config_ID": c.config_id,
                "Transcript": c.transcript_id,
                "TSS": c.canonical_tss,
                "MANE": "MANE Select" if c.is_mane_select else ("MANE Plus Clinical" if c.is_mane_plus_clinical else "No"),
                "Window": c.window_name,
                "Promoter_Range": f"[{c.promoter_start}, {c.promoter_end})",
                "Experiment": c.encode_experiment_accession,
                "File": c.encode_file_accession,
                "Biosample": c.biosample,
                "Treatment": c.treatment_context,
                "Modification": c.genetic_modification_context,
                "Tier": c.processing_tier,
                "Auth_Status": c.auth_status,
                "Overlap_bp": c.overlap_bp,
                "Overlap_pct": f"{c.overlap_percentage:.1f}%",
                "Evidence_Status": status_val,
                "Exclusion_Reason": c.exclusion_reason or ""
            })
        return pd.DataFrame(rows)

    def to_markdown(self) -> str:
        """Renders the matrix and causal disagreement diagnoses as GitHub Flavored Markdown."""
        stats = self.get_summary_statistics()
        analysis = self.analyze()

        lines = [
            f"# Regulatory Evidence Disagreement Matrix: {self.gene_symbol} / {self.tf_target}",
            "",
            "### Evaluation Summary",
            f"- **Target Locus**: `{self.gene_symbol}` | **TF Target**: `{self.tf_target}`",
            f"- **Total Configurations**: {stats['total_configurations']}",
            f"- **Independent Experiments**: {stats['unique_experiments_count']}",
            f"- **Evaluated Transcripts**: {stats['unique_transcripts_count']} | **Promoter Windows**: {stats['unique_windows_count']}",
            f"- **Evidence Breakdown**: {stats['positive_configurations_count']} Verified Positive, {stats['negative_configurations_count']} Verified No Peak, {stats['excluded_configurations_count']} Excluded",
            f"- **Overall Concordance Status**: `{stats['overall_status']}`",
            "",
            "> [!IMPORTANT]",
            f"> **Scientific Interpretation**: {stats['scientific_interpretation']}",
            "",
            "### Configuration Evidence Grid",
            "| Config ID | Transcript | TSS | Window | Range | Experiment | File | Context | Auth Status | Overlap | Evidence |",
            "|:---|:---|:---:|:---|:---:|:---|:---|:---|:---:|:---:|:---:|"
        ]

        for c in self.configurations:
            ctx_summary = f"{c.biosample} ({c.treatment_context})"
            range_str = f"{c.promoter_start}-{c.promoter_end}"
            if c.evidence_status == EvidenceInterpretation.VERIFIED_POSITIVE_OBSERVATION:
                status_badge = "🟢 VERIFIED POSITIVE"
            elif c.evidence_status == EvidenceInterpretation.VERIFIED_NO_PEAK_DETECTED:
                status_badge = "⚪ VERIFIED NO PEAK"
            else:
                status_badge = "🔴 EXCLUDED / UNVERIFIED"

            overlap_str = f"{c.overlap_bp} bp ({c.overlap_percentage:.1f}%)" if c.has_overlap else "0 bp"
            lines.append(
                f"| `{c.config_id}` | `{c.transcript_id}` | {c.canonical_tss:,} | {c.window_name} | "
                f"`{range_str}` | `{c.encode_experiment_accession}` | `{c.encode_file_accession}` | "
                f"{ctx_summary} | `{c.auth_status}` | {overlap_str} | **{status_badge}** |"
            )

        if analysis["diagnoses"]:
            lines.extend([
                "",
                "### Detailed Cause-of-Disagreement Diagnoses",
                "| Pair | Primary Causes | Smallest Observed Difference | Scientific Explanation | Attribution Caveat |",
                "|:---|:---|:---|:---|:---|"
            ])
            for d in analysis["diagnoses"]:
                causes_badge = "<br>".join([f"`{c}`" for c in d["primary_causes"]])
                pair_lbl = f"`{d['config_a_id']}` vs `{d['config_b_id']}`"
                caveat = d["causal_attribution_caveat"] or "N/A"
                lines.append(
                    f"| {pair_lbl} | {causes_badge} | {d['smallest_observed_difference']} | "
                    f"{d['scientific_explanation']} | {caveat} |"
                )

        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        """Serializes the matrix and evaluation summary to a dictionary."""
        stats = self.get_summary_statistics()
        analysis = self.analyze()
        return {
            "gene_symbol": self.gene_symbol,
            "tf_target": self.tf_target,
            "evaluation_summary": {
                "total_configurations": stats["total_configurations"],
                "positive_count": stats["positive_configurations_count"],
                "negative_count": stats["negative_configurations_count"],
                "excluded_count": stats["excluded_configurations_count"],
                "overall_status": stats["overall_status"],
                "scientific_interpretation": stats["scientific_interpretation"]
            },
            "configurations": [c.to_dict() for c in self.configurations],
            "analysis": analysis
        }

    def to_json(self, indent: int = 2) -> str:
        """Serializes the matrix to a JSON string."""
        return json.dumps(self.to_dict(), indent=indent, default=str)

