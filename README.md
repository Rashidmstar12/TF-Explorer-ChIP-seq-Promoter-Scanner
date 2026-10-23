# TF-Explorer: Epigenetic Promoter Scanner & Regulatory Disagreement Engine

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Version 3.1.0rc1](https://img.shields.io/badge/version-3.1.0rc1-green.svg)](https://github.com/Rashidmstar12/TF-Explorer-ChIP-seq-Promoter-Scanner)

**TF-Explorer** is an open-source computational genomics toolkit and interactive analysis suite designed to scan, cross-validate, and triage transcription factor (TF) binding evidence across human gene promoters.

By integrating authoritative **ENCODE ChIP-seq** peak datasets, **JASPAR PSSM** motif scans, and **MANE Select v1.3** transcript annotations, TF-Explorer addresses the pervasive challenge of discordant regulatory evidence in genomic studies.

---

## 🧬 Scientific Problem & The Disagreement Engine

Genomic researchers frequently encounter conflicting regulatory evidence:
- An ENCODE ChIP-seq peak is detected in a promoter, but no canonical motif sequence exists.
- A high-scoring motif exists in the DNA sequence, but no ChIP-seq signal is observed.
- Two studies report opposing binding results because they employed different promoter window definitions (e.g., -1000 bp vs. -2000 bp from TSS) or alternative transcript isoforms.

TF-Explorer addresses this challenge with the **Regulatory Evidence Disagreement Engine**, which automates evidence triage and discordance classification across a **Seven-Category Regulatory Discordance Taxonomy**:

1. **`DIRECT_COINCIDENT_OCCUPANCY`**: Coincident genomic evidence: coincident ChIP-seq peak and JASPAR motif observed.
2. **`MOTIF_UNSUPPORTED_CHIP_PEAK`**: ChIP-seq peak observed without coincident canonical sequence motif (consistent with indirect tethering or co-factor recruitment).
3. **`UNBOUND_CANONICAL_MOTIF`**: Canonical motif present without overlapping peak (consistent with closed chromatin, DNA methylation, or absence of TF expression).
4. **`PROMOTER_WINDOW_DEPENDENT`**: Discrepancy caused solely by divergent promoter window arithmetic relative to the TSS.
5. **`TRANSCRIPT_ISOFORM_SPECIFIC`**: Discrepancy arising from alternative transcription start sites (TSS) between transcript isoforms.
6. **`CELL_TYPE_SPECIFIC_BINDING`**: Differential occupancy observed across distinct cellular contexts (e.g., K562 vs. HepG2).
7. **`TECHNICAL_ARTIFACT_OR_QUALITY_DISCORDANCE`**: Discrepancy driven by sequencing depth, peak calling thresholds ($p$-value / $q$-value), or processing pipelines.

---

## 💻 Installation

### Option 1: Install from Pre-Built Wheel (Recommended)

```bash
pip install tf_explorer-3.1.0rc1-py3-none-any.whl
```

### Option 2: Install from Source

```bash
git clone https://github.com/Rashidmstar12/TF-Explorer-ChIP-seq-Promoter-Scanner.git
cd TF-Explorer-ChIP-seq-Promoter-Scanner
pip install .
```

For development and test dependencies:
```bash
pip install -e ".[dev,test]"
```

---

## 🚀 Quick Start

### 1. Minimal Offline Demonstration (<0.1s Execution)

Verify core coordinate and dispute triage logic immediately without internet access:

```bash
python examples/quickstart_promoter.py
python examples/quickstart_disagreement.py
```

### 2. Command-Line Interface (CLI)

Scan ENCODE ChIP-seq peaks and JASPAR motifs across any human gene promoter:

```bash
tf-explorer \
  --gene TP53 \
  --tf-list "E2F1,CTCF" \
  --jaspar-ids "MA0024.1,MA0139.1" \
  --promoter-up 2000 \
  --promoter-down 500 \
  --out results/
```

#### Key CLI Arguments:
- `--gene`: Target gene symbol (e.g., `TP53`, `CDK4`, `GAPDH`).
- `--tf-list`: Comma-separated list of transcription factors to scan.
- `--jaspar-ids`: Comma-separated JASPAR matrix IDs for PWM scoring.
- `--genome`: Reference genome assembly (default: `hg38` / `GRCh38`).
- `--promoter-up`: Base pairs upstream of TSS (default: `2000`).
- `--promoter-down`: Base pairs downstream of TSS (default: `500`).
- `--threshold`: JASPAR PWM log-odds score threshold (default: `8.0`).
- `--out`: Output directory for generated CSV tables and visualizations.

### 3. Python API

Programmatically classify regulatory discordance in custom analysis pipelines:

```python
from tf_explorer import calculate_promoter_coordinates, DisagreementEngine, EvidenceConfiguration, EvidenceInterpretation

# 1. Expand strand-aware promoter boundaries
start, end = calculate_promoter_coordinates(tss=7676594, strand="-", upstream=2000, downstream=500)
print(f"TP53 promoter interval: chr17:{start}-{end}")

# 2. Triage conflicting regulatory evidence configurations
engine = DisagreementEngine()
# evaluate_pair diagnoses causes across 7 discordance categories
```

### 4. Interactive Streamlit Web Interface

Launch the full interactive visualization suite:

```bash
streamlit run tf_explorer/app.py
```

---

## 🔬 Methodological Boundaries & Evidence Triage Caveats

- **Evidence Triage, Not Physical Assay:** TF-Explorer is a computational evidence-triage and discordance classification tool. Inferred binding states and dispute taxonomy diagnoses represent computational hypotheses based on public genomic data and user-specified geometric criteria, not direct physical verification of in vivo binding. Empirical experimental validation (e.g., EMSA, CUT&RUN, or functional reporter assays) remains essential.
- **Reference Genome:** GRCh38 / hg38. All coordinate calculations enforce 0-based half-open interval conventions $[start, end)$.
- **Strand Geometry:** On the negative strand ($-$), promoter windows expand to numerically higher coordinates $[TSS - down, TSS + up]$ ensuring $start < end$ globally.
- **Transcript Canonical Authority:** Transcript TSS coordinates are verified against MANE Select v1.3 (Ensembl / NCBI consensus).
- **Public REST Clients:** External queries (ENCODE, Ensembl, UCSC, GTEx, STRING) use standard pure-Python HTTP clients with enforced timeouts and retry backoffs. No private shell tools or local agent skills are required.

---

## 🧪 Running the Test Suite

Execute the automated test suite with pytest:

```bash
pytest tests/ -v
```

---

## 📄 License & Governance

TF-Explorer is licensed under the **OSI-Approved MIT License**. See [LICENSE](LICENSE) for details.

Contributions are welcome! Please review our [Contributing Guidelines](CONTRIBUTING.md) and [Code of Conduct](CODE_OF_CONDUCT.md).

For inquiries or bug reports, open an issue on GitHub or contact `rashidmstar@pondiuni.ac.in`.

---

## 📊 Research Impact & JOSS Eligibility Status

**Current JOSS Research Status:** JOSS RESEARCH IMPACT: NOT YET DEMONSTRATED

- **Computational Tooling vs. Empirical Discovery:** TF-Explorer v3.1.0rc1 provides computational triage algorithms, geometric coordinate reconciliation, and authentication contracts. A passing automated test suite verifies expected code behaviors under specified conditions, not the absence of every conceivable defect or physical biological binding.
- **Protected Prospective Benchmark:** Prospective biological benchmark cases remain strictly blinded and frozen to protect benchmark integrity and prevent post-hoc overfitting.
- **Synthetic Demonstrations:** Quickstart scripts and test dispute scenarios represent synthetic demonstrations for testing and tutorial purposes, not empirical biological findings.
- **Future JOSS Submission:** Formal submission to the Journal of Open Source Software (JOSS) is deferred until after public release, community review, independent human-led biological research applications, and complete manuscript-codebase synchronization.

