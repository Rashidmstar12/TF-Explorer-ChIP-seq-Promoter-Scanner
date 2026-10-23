"""
Conclusion Robustness and Sensitivity Analyzer for TF-Explorer.
Phase 3B Pilot Implementation.

Evaluates whether promoter-associated regulatory observations remain stable or fragile
under alternative transcript models, promoter window definitions, and experimental datasets.
Does NOT generate artificial black-box confidence scores.
"""

import os
import logging
from enum import Enum
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field, asdict

logger = logging.getLogger(__name__)


class RobustnessClassification(str, Enum):
    STABLE_WITHIN_EVALUATED_CHOICES = "STABLE_WITHIN_EVALUATED_CHOICES"
    SENSITIVE_TO_TRANSCRIPT = "SENSITIVE_TO_TRANSCRIPT"
    SENSITIVE_TO_EXPERIMENT = "SENSITIVE_TO_EXPERIMENT"
    SENSITIVE_TO_PARAMETERS = "SENSITIVE_TO_PARAMETERS"
    CONTEXT_DEPENDENT = "CONTEXT_DEPENDENT"
    DISCORDANT_EVIDENCE = "DISCORDANT_EVIDENCE"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass
class ConfigurationEvaluation:
    """Evaluation result for one specific analytical configuration."""
    transcript_id: str
    transcript_name: str
    is_canonical: bool
    tss: int
    strand: str
    window_upstream: int
    window_downstream: int
    promoter_start: int
    promoter_end: int
    experiment_id: str
    biosample: str
    quality_tier: str
    has_peak_overlap: bool
    overlapping_peak_count: int
    overlapping_peaks: List[Dict[str, Any]] = field(default_factory=list)
    has_motif_match: Optional[bool] = None
    motif_match_count: int = 0
    motif_matches: List[Dict[str, Any]] = field(default_factory=list)
    data_status: str = "DATA_AVAILABLE"
    biological_context: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class RobustnessReport:
    """Summary of multiverse parameter sensitivity analysis."""
    gene_symbol: str
    target_tf: str
    biosample: str
    classification: RobustnessClassification
    support_ratio: float
    total_valid_configurations: int
    supporting_configurations: int
    experiment_support_fraction: str
    transcript_support: Dict[str, bool]
    window_support: Dict[str, bool]
    configurations: List[ConfigurationEvaluation]
    operational_rationale: str
    negative_inference_guardrail: str
    biological_contexts: Dict[str, Any] = field(default_factory=dict)
    context_warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["classification"] = self.classification.value
        d["configurations"] = [c.to_dict() for c in self.configurations]
        return d

    def to_researcher_markdown(self) -> str:
        """Renders an intuitive, scientifically complete Markdown evidence report."""
        lines = []
        lines.append(f"# Regulatory Evidence Robustness Dossier: {self.gene_symbol} — {self.target_tf}")
        lines.append(f"**Target Gene**: `{self.gene_symbol}` | **Transcription Factor**: `{self.target_tf}` | **Biosample**: `{self.biosample}`")
        lines.append(f"**Multiverse Robustness Classification**: `{self.classification.value}`")
        lines.append(f"**Parameter Grid Stability**: `{self.supporting_configurations}/{self.total_valid_configurations}` ({self.support_ratio:.1%}) | **Independent Experiments**: `{self.experiment_support_fraction}`")
        lines.append("")
        lines.append(f"> [!NOTE]\n> **{self.negative_inference_guardrail}**")
        lines.append("")

        if self.context_warnings:
            lines.append("### Experimental Context & Perturbation Alerts")
            for w in self.context_warnings:
                lines.append(f"- ⚠️ **Alert**: {w}")
            lines.append("")

        if self.biological_contexts:
            lines.append("### Evaluated Experiments & Contexts")
            lines.append("| Experiment Accession | Biosample | Context Summary | Treatments | Modifications / Tags |")
            lines.append("| :--- | :--- | :--- | :--- | :--- |")
            for exp_id, ctx in self.biological_contexts.items():
                bio = ctx.get("biosample_term", self.biosample)
                summ = ctx.get("biosample_summary", "Standard culture")
                tr = ", ".join([t.get("treatment_term_name", str(t)) if isinstance(t, dict) else str(t) for t in ctx.get("treatments", [])]) or "None (Untreated)"
                tags = []
                for gm in ctx.get("genetic_modifications", []):
                    if isinstance(gm, dict):
                        for tag in gm.get("tags", []):
                            tags.append(tag.get("name", "tag"))
                tag_str = ", ".join(tags) or "None (Wild-type)"
                lines.append(f"| `{exp_id}` | {bio} | {summ} | {tr} | {tag_str} |")
            lines.append("")

        lines.append("### Multiverse Parameter Sensitivity Grid")
        lines.append("| Transcript | TSS | Strand | Window | Promoter Interval | Overlapping Peaks | Overlap bp | Status |")
        lines.append("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |")
        for c in self.configurations:
            canon = " (Canonical)" if c.is_canonical else ""
            pk_str = f"{c.overlapping_peak_count} peak(s)" if c.has_peak_overlap else "0 peaks"
            ov_bp_str = f"{c.overlapping_peaks[0]['overlap_bp']} bp" if c.overlapping_peaks else "0 bp"
            lines.append(
                f"| `{c.transcript_id}`{canon} | {c.tss:,} | `{c.strand}` | `-{c.window_upstream}/+{c.window_downstream}` | "
                f"`{c.promoter_start:,}–{c.promoter_end:,}` | {pk_str} | {ov_bp_str} | `{c.data_status}` |"
            )
        lines.append("")

        lines.append("### Operational Rationale")
        lines.append(self.operational_rationale)
        lines.append("")
        return "\n".join(lines)



class RobustnessAnalyzer:
    """
    Multiverse sensitivity analyzer evaluating promoter regulatory evidence
    across transcripts, promoter window sizes, and independent experiments.
    """

    def __init__(self):
        pass

    @staticmethod
    def calculate_promoter_coordinates(tss: int, strand: str, up_bp: int, down_bp: int) -> Tuple[int, int]:
        """
        Calculates 1-based closed genomic coordinates for a promoter window.
        Correctly handles strand orientation:
        - For (+) strand: 5' upstream is < TSS, downstream is > TSS: [TSS - up, TSS + down]
        - For (-) strand: 5' upstream is > TSS, downstream is < TSS: [TSS - down, TSS + up]
        """
        if strand == "+":
            return tss - up_bp, tss + down_bp
        else:
            return tss - down_bp, tss + up_bp

    @staticmethod
    def check_peak_overlap(peak_start: int, peak_end: int, prom_start: int, prom_end: int) -> Tuple[bool, int]:
        """
        Calculates overlap between peak interval and promoter interval.
        Returns (is_overlapping, overlap_bp).
        """
        overlap_start = max(peak_start, prom_start)
        overlap_end = min(peak_end, prom_end)
        if overlap_start < overlap_end:
            return True, overlap_end - overlap_start
        return False, 0

    def evaluate(
        self,
        gene_symbol: str,
        target_tf: Optional[str] = None,
        biosample: str = "Unknown",
        transcripts: Optional[List[Dict[str, Any]]] = None,
        experiment_peak_sets: Optional[List[Dict[str, Any]]] = None,
        promoter_windows: Optional[List[Tuple[int, int]]] = None,
        motif_scanner_fn: Optional[Any] = None,
        biological_contexts: Optional[Dict[str, Any]] = None,
        tf: Optional[str] = None,
        windows: Optional[List[Tuple[int, int]]] = None
    ) -> RobustnessReport:
        """
        Evaluates regulatory evidence across a grid of:
        Transcripts x Promoter Windows x Experiments.

        Parameters:
        - gene_symbol: Gene of interest.
        - target_tf: Transcription factor of interest.
        - biosample: Cellular context.
        - transcripts: List of transcripts from Ensembl with keys:
                       'id', 'name', 'tss', 'strand', 'is_canonical'.
        - experiment_peak_sets: List of dicts representing independent experiments:
                                {'experiment_id': str, 'quality_tier': str, 'peaks': List[Dict]}.
        - promoter_windows: List of (up_bp, down_bp) tuples, e.g. [(500, 100), (1000, 200), (2000, 500)].
        - motif_scanner_fn: Optional callable (chrom, start, end, strand) -> List[Dict].
        - biological_contexts: Optional dict of experiment biological contexts.
        """
        if transcripts is None:
            transcripts = []
        if experiment_peak_sets is None:
            experiment_peak_sets = []
        actual_tf = target_tf or tf or "Unknown_TF"
        if promoter_windows is None:
            promoter_windows = windows or [(500, 100), (1000, 200), (2000, 500)]

        guardrail = (
            "STATISTICAL GUARDRAIL: The support ratio reflects geometric stability across tested parameter choices, "
            "NOT an in vivo biological probability of transcriptional regulation. A support ratio of 1.0 does not prove "
            "activation, and discordance between experiments may reflect biological cell-state or antibody differences."
        )

        if not transcripts:
            return RobustnessReport(
                gene_symbol=gene_symbol,
                target_tf=target_tf,
                biosample=biosample,
                classification=RobustnessClassification.INSUFFICIENT_DATA,
                support_ratio=0.0,
                total_valid_configurations=0,
                supporting_configurations=0,
                experiment_support_fraction="0/0",
                transcript_support={},
                window_support={},
                configurations=[],
                operational_rationale="No annotated transcript models available for gene.",
                negative_inference_guardrail=guardrail
            )

        if not experiment_peak_sets:
            return RobustnessReport(
                gene_symbol=gene_symbol,
                target_tf=target_tf,
                biosample=biosample,
                classification=RobustnessClassification.INSUFFICIENT_DATA,
                support_ratio=0.0,
                total_valid_configurations=0,
                supporting_configurations=0,
                experiment_support_fraction="0/0",
                transcript_support={},
                window_support={},
                configurations=[],
                operational_rationale=f"No experimental ChIP-seq datasets available for {target_tf} in {biosample}.",
                negative_inference_guardrail=guardrail
            )

        config_evals: List[ConfigurationEvaluation] = []
        transcript_support: Dict[str, List[bool]] = {t['id']: [] for t in transcripts}
        window_support_map: Dict[str, List[bool]] = {f"-{w[0]}/+{w[1]}": [] for w in promoter_windows}
        exp_support_map: Dict[str, List[bool]] = {e['experiment_id']: [] for e in experiment_peak_sets}

        for t in transcripts:
            t_id = t['id']
            t_name = t.get('name', t_id)
            is_canon = t.get('is_canonical', False)
            tss = t['tss']
            strand = t['strand']

            for w_up, w_down in promoter_windows:
                w_label = f"-{w_up}/+{w_down}"
                p_start, p_end = self.calculate_promoter_coordinates(tss, strand, w_up, w_down)

                for exp_entry in experiment_peak_sets:
                    exp_id = exp_entry['experiment_id']
                    q_tier = exp_entry.get('quality_tier', 'optimal')
                    peaks = exp_entry.get('peaks', [])

                    overlapping_peaks = []
                    t_chrom = t.get('chrom') or t.get('chromosome')
                    for pk in peaks:
                        pk_chrom = pk.get('chrom') or pk.get('chromosome')
                        if pk_chrom and t_chrom:
                            norm_p = str(pk_chrom) if str(pk_chrom).startswith('chr') else f"chr{pk_chrom}"
                            norm_t = str(t_chrom) if str(t_chrom).startswith('chr') else f"chr{t_chrom}"
                            if norm_p != norm_t:
                                continue
                        pk_start = pk.get('start', 0)
                        pk_end = pk.get('end', 0)
                        is_ov, ov_bp = self.check_peak_overlap(pk_start, pk_end, p_start, p_end)
                        if is_ov:
                            overlapping_peaks.append({
                                'peak_id': pk.get('id', f"peak_{pk_start}_{pk_end}"),
                                'start': pk_start,
                                'end': pk_end,
                                'overlap_bp': ov_bp,
                                'signal_value': pk.get('signal_value', 0.0),
                                'q_value': pk.get('q_value', 0.0)
                            })

                    has_peak = len(overlapping_peaks) > 0

                    motif_matches = []
                    has_motif = None
                    if motif_scanner_fn is not None and has_peak:
                        try:
                            motif_matches = motif_scanner_fn(t.get('chrom', ''), p_start, p_end, strand)
                            has_motif = len(motif_matches) > 0
                        except Exception as e_m:
                            logger.warning(f"Motif scan failed for {t_id} {w_label}: {e_m}")

                    ce = ConfigurationEvaluation(
                        transcript_id=t_id,
                        transcript_name=t_name,
                        is_canonical=is_canon,
                        tss=tss,
                        strand=strand,
                        window_upstream=w_up,
                        window_downstream=w_down,
                        promoter_start=p_start,
                        promoter_end=p_end,
                        experiment_id=exp_id,
                        biosample=biosample,
                        quality_tier=q_tier,
                        has_peak_overlap=has_peak,
                        overlapping_peak_count=len(overlapping_peaks),
                        overlapping_peaks=overlapping_peaks,
                        has_motif_match=has_motif,
                        motif_match_count=len(motif_matches),
                        motif_matches=motif_matches,
                        data_status="PEAK_DETECTED" if has_peak else "ASSAYED_NO_PEAK_DETECTED",
                        biological_context=exp_entry.get("biological_context", {})
                    )
                    config_evals.append(ce)

                    transcript_support[t_id].append(has_peak)
                    window_support_map[w_label].append(has_peak)
                    exp_support_map[exp_id].append(has_peak)

        total_configs = len(config_evals)
        supporting_configs = sum(1 for c in config_evals if c.has_peak_overlap)
        support_ratio = supporting_configs / total_configs if total_configs > 0 else 0.0

        # Summarize by dimension
        final_tx_support = {tid: (any(vals) if vals else False) for tid, vals in transcript_support.items()}
        final_w_support = {w: (any(vals) if vals else False) for w, vals in window_support_map.items()}
        
        # Experiment-level agreement grouped by unique biological experiment ID (do not conflate reprocessing files)
        unique_exp_ids = list(dict.fromkeys(e['experiment_id'] for e in experiment_peak_sets))
        total_unique_exps = len(unique_exp_ids)
        exp_supported_count = sum(1 for exp_id in unique_exp_ids if any(exp_support_map.get(exp_id, [])))
        exp_fraction_str = f"{exp_supported_count}/{total_unique_exps}"

        # Collect biological contexts and warnings
        exp_contexts = {e['experiment_id']: e.get('biological_context', {}) for e in experiment_peak_sets}
        if biological_contexts:
            exp_contexts.update(biological_contexts)
        all_warnings = []
        has_perturbed_context = False
        for exp_id, ctx in exp_contexts.items():
            if ctx.get("is_treated") or ctx.get("is_genetically_modified") or ctx.get("treatments") or ctx.get("genetic_modifications"):
                has_perturbed_context = True
            warnings_list = ctx.get('context_warnings', []) or ctx.get('notes', [])
            for w in warnings_list:
                if w not in all_warnings:
                    all_warnings.append(w)

        # Determine Operational Classification
        classification = RobustnessClassification.STABLE_WITHIN_EVALUATED_CHOICES
        rationale_parts = []

        # Check Discordance between experiments
        if total_unique_exps > 1 and exp_supported_count > 0 and exp_supported_count < total_unique_exps:
            # Some experiments found peak, others did not
            classification = RobustnessClassification.SENSITIVE_TO_EXPERIMENT
            rationale_parts.append(
                f"Peak occupancy is sensitive to experiment choice: supported in {exp_supported_count} of {total_unique_exps} independent experiments."
            )
        elif supporting_configs == 0:
            classification = RobustnessClassification.STABLE_WITHIN_EVALUATED_CHOICES
            rationale_parts.append("Consistent negative binding: no peak overlaps observed across any evaluated transcript or window.")
        elif all(c.has_peak_overlap for c in config_evals):
            classification = RobustnessClassification.STABLE_WITHIN_EVALUATED_CHOICES
            rationale_parts.append("Observation is stable: peak overlap observed across all evaluated transcripts, windows, and experiments.")
            if has_perturbed_context and all_warnings:
                rationale_parts.append("Note: While parameter stability is maintained across evaluated choices, the observation reflects an experimental context with biological limitations.")
        else:
            # Variation exists. Check if variation is strictly driven by transcripts
            tx_has_diff = len(set(final_tx_support.values())) > 1
            w_has_diff = any(not all(window_support_map[w]) for w in window_support_map)

            if tx_has_diff:
                classification = RobustnessClassification.SENSITIVE_TO_TRANSCRIPT
                canon_tx = [t['id'] for t in transcripts if t.get('is_canonical')]
                canon_status = final_tx_support.get(canon_tx[0], False) if canon_tx else False
                rationale_parts.append(
                    f"Observation depends on transcript choice. Supported in {sum(final_tx_support.values())} of {len(transcripts)} transcripts "
                    f"(Canonical transcript {canon_tx}: {'SUPPORTED' if canon_status else 'UNSUPPORTED'})."
                )
                # Granular geometric explanation of transcript differences
                sup_txs = [t for t in transcripts if final_tx_support.get(t['id'])]
                unsup_txs = [t for t in transcripts if not final_tx_support.get(t['id'])]
                if sup_txs and unsup_txs:
                    s_t = sup_txs[0]
                    u_t = unsup_txs[0]
                    disp = abs(s_t['tss'] - u_t['tss'])
                    rationale_parts.append(
                        f"Peak overlap is observed for transcript {s_t['id']} (TSS {s_t['tss']:,}) but not transcript {u_t['id']} (TSS {u_t['tss']:,}) "
                        f"because transcription start sites differ by {disp:,} bp, placing the peak outside {u_t['id']}'s evaluated promoter."
                    )
            elif w_has_diff:
                classification = RobustnessClassification.SENSITIVE_TO_PARAMETERS
                rationale_parts.append(
                    "Observation depends on promoter window definition: peak overlap is present in expanded distal promoter windows but absent in core proximal windows."
                )
            else:
                classification = RobustnessClassification.SENSITIVE_TO_PARAMETERS
                rationale_parts.append("Observation varies across parameter configurations.")

        return RobustnessReport(
            gene_symbol=gene_symbol,
            target_tf=actual_tf,
            biosample=biosample,
            classification=classification,
            support_ratio=round(support_ratio, 4),
            total_valid_configurations=total_configs,
            supporting_configurations=supporting_configs,
            experiment_support_fraction=exp_fraction_str,
            transcript_support=final_tx_support,
            window_support=final_w_support,
            configurations=config_evals,
            operational_rationale=" ".join(rationale_parts),
            negative_inference_guardrail=guardrail,
            biological_contexts=exp_contexts,
            context_warnings=all_warnings
        )

    analyze = evaluate

