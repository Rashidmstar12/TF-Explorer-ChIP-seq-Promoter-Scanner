# Release Source Identity Correction Report: Forensic ENCODE Accession Resolution

**Document ID:** `RELEASE-SRC-IDENTITY-001`  
**Date:** 2026-10-23  
**Status:** **DEFECT CONFIRMED, DOCUMENTATION RETRACTED, AUTOMATED CHECK ENFORCED**  
**Component:** `tf_explorer/integrity.py`, `tf_explorer/encode_client.py`, `FINAL_RESEARCH_IMPACT_STATUS.md`

---

## 1. Executive Summary & Defect Declaration

In Phase 4C-8R, `FINAL_RESEARCH_IMPACT_STATUS.md` proposed a future real-data empirical research plan to demonstrate research utility for JOSS:
- **Claimed Target:** Estrogen Receptor (*ESR1*) in MCF-7 cells.
- **Claimed Locus:** Trefoil Factor 1 (*TFF1*) promoter.
- **Asserted Experiment Accession:** `ENCSR000BST`
- **Asserted File Accession:** `ENCFF001UQF`

An independent audit of authoritative ENCODE metadata established two critical source attribution errors:
1. **`ENCSR000BST` does NOT target ESR1.** The authoritative ENCODE portal record confirms that `ENCSR000BST` is a ChIP-seq experiment targeting **GATA3** in MCF-7 cells.
2. **`ENCFF001UQF` does NOT belong to `ENCSR000BST`.** Authoritative file metadata confirms that `ENCFF001UQF` belongs to dataset `/experiments/ENCSR000BQF/`, which targets **TAF1** in **SK-N-SH** cells (neuroblastoma) on assembly **hg19**.

Furthermore, prior project documentation in Phase 4C-8 had erroneously linked `ENCFF001UQF` with SPI1 in K562 cells.

These errors demonstrate a recurring defect: **manual assertion of accessions in documentation without programmatic verification against authoritative metadata.**

This report documents the forensic reconciliation, traces every occurrence of incorrect accessions, defines an automated consistency check, outlines the mandatory 6-gate verification protocol, and formally declares the proposed study **NOT AUTHORIZED FOR EXECUTION** until verified.

---

## 2. Authoritative ENCODE Metadata Verification

### 2.1 Experiment `ENCSR000BST`
* **Source:** `https://www.encodeproject.org/experiments/ENCSR000BST/`
* **Query Timestamp:** 2026-10-23
* **Official Attributes:**
  - **Accession:** `ENCSR000BST`
  - **Target Label:** `GATA3` (GATA binding protein 3, human)
  - **Target Gene:** `GATA3` (NCBI Gene ID: 2625, Ensembl: ENSG00000107485)
  - **Biosample:** `MCF-7` (human breast adenocarcinoma)
  - **Assay Title:** `TF ChIP-seq`
  - **Description:** `"GATA3 ChIP-seq protocol v042211.1 on human MCF-7"`
  - **Files in Dataset:** 32 files (all targeting GATA3)
  - **Presence of `ENCFF001UQF`:** **ABSENT (`False`)**
* **Finding:** Attributing `ENCSR000BST` to ESR1 is a scientific error.

### 2.2 File `ENCFF001UQF`
* **Source:** `https://www.encodeproject.org/files/ENCFF001UQF/`
* **Query Timestamp:** 2026-10-23
* **Official Attributes:**
  - **File Accession:** `ENCFF001UQF`
  - **Dataset / Experiment:** `/experiments/ENCSR000BQF/`
  - **File Format:** `bed` (`narrowPeak`)
  - **Output Type:** `peaks`
  - **Genome Assembly:** `hg19` (GRCh37)
  - **File Status:** `released`
  - **MD5 Checksum:** `5e76dd10fb62e93d8e58c6252694b757`
* **Parent Dataset Attributes (`ENCSR000BQF`):**
  - **Experiment Accession:** `ENCSR000BQF`
  - **Target Label:** `TAF1` (TATA-box binding protein associated factor 1)
  - **Biosample:** `SK-N-SH` (human neuroblastoma)
  - **Description:** `"TAF1 ChIP-seq protocol v041610.1 on human SK-N-SH"`
* **Finding:** `ENCFF001UQF` targets TAF1 in SK-N-SH on hg19. It does not belong to `ENCSR000BST`, does not target ESR1, does not target GATA3, and does not originate from MCF-7 or K562.

---

## 3. Search and Audit of Reused Accessions Across the Codebase

A search across the entire project repository identified all occurrences of `ENCSR000BST` and `ENCFF001UQF`:

| File / Document | Accession Cited | Context / Prior Claim | Audit Finding |
| :--- | :--- | :--- | :--- |
| `FINAL_RESEARCH_IMPACT_STATUS.md` | `ENCSR000BST` | Claimed as ESR1 ChIP-seq in MCF-7 | **ERRONEOUS:** Retracted. Targets GATA3. |
| `FINAL_RESEARCH_IMPACT_STATUS.md` | `ENCFF001UQF` | Claimed as ESR1 peak file in MCF-7 | **ERRONEOUS:** Retracted. Belongs to TAF1 / SK-N-SH (`ENCSR000BQF`). |
| `FINAL_RELEASE_ALLOWLIST_VERIFICATION.md` | `ENCFF001UQF` | Listed as hypothetical subset in `reference_data/` | **FALSE STATEMENT:** File never existed in release package. Retracted. |
| `FINAL_AUTHOR_RELEASE_APPROVAL.md` | `ENCFF001UQF` | Cited in author research validation plan | **ERRONEOUS:** Retracted. |
| `FINAL_JOSS_ELIGIBILITY_MATRIX.md` | `ENCFF001UQF` | Cited in roadmap unblocking table | **ERRONEOUS:** Retracted. |
| `PHASE4C8_SOURCE_TRUTH_CORRECTION.md` | `ENCFF001UQF` | Associated with K562 / SPI1 context | **ERRONEOUS:** Retracted. |

**Code Integrity Confirmation:** Zero occurrences of `ENCSR000BST` or `ENCFF001UQF` exist within the core library source (`tf_explorer/`) or the 18 public tests (`tests/`). Both accessions were confined exclusively to audit reports, proposals, and narrative text.

---

## 4. Root Cause Analysis: Why Biological Attribution Errors Recurred

1. **Reliance on Informal Recollection:** In drafting research proposals in Phase 4C-7 and 4C-8R, accession identifiers were cited based on informal recollection or draft notes without running the code's verification engine.
2. **Decoupled Documentation Authoring:** Documentation authoring was not gated by the programmatic validation pipeline that guards the core software.
3. **Cross-Contamination of Identifiers:** When GATA1 (`ENCSR000BGW`) was corrected in Phase 4C-8, the author substituted another breast cancer dataset from memory without verifying that `ENCSR000BST` was GATA3 rather than ESR1.

---

## 5. Mandatory 6-Gate Source Verification Protocol

To prevent recurrence, no biological accession may enter any research protocol, publication report, or data ingestion pipeline without satisfying all six gates:

```
+-----------------------------------------------------------------------------------+
|                           MANDATORY 6-GATE PROTOCOL                               |
+===================================================================================+
| Gate 1: Target Verification      | target.label == expected_target                |
| Gate 2: Biosample Verification   | biosample_ontology.term_name == expected_cell  |
| Gate 3: File Membership          | file.dataset == experiment.accession           |
| Gate 4: Assembly Verification    | file.assembly == expected_assembly (GRCh38)    |
| Gate 5: Output Type Verification | file.output_type in accepted_peak_types        |
| Gate 6: Checksum Verification    | file.md5sum == primary_file_md5                |
+-----------------------------------------------------------------------------------+
```

---

## 6. Automated Source-Identity Consistency Check

The codebase already contains an automated 12-point integrity engine: [`EvidenceIntegrityEngine`](file:///c:/Users/rashi/Documents/TF%20EXPLORER/tf_explorer/integrity.py).

When evaluated against the mismatched pair (`ENCFF001UQF`, `ENCSR000BST`), `EvidenceIntegrityEngine.verify()` immediately detects and rejects the record fail-closed:

```python
from tf_explorer.integrity import EvidenceIntegrityEngine

engine = EvidenceIntegrityEngine()
res = engine.verify(
    file_accession="ENCFF001UQF",
    experiment_accession="ENCSR000BST",
    expected_tf="ESR1",
    expected_biosample="MCF-7",
    expected_assembly="GRCh38"
)

# Output:
# res.is_valid == False
# res.status == ValidationStatusCode.FILE_EXPERIMENT_MISMATCH
# res.failure_code == "FILE_EXPERIMENT_MISMATCH"
# res.message == "File 'ENCFF001UQF' belongs to parent experiment 'ENCSR000BQF', not 'ENCSR000BST'."
```

If `ENCSR000BST` is checked against requested TF `ESR1`, the engine detects:
```python
# res.failure_code == "TF_TARGET_MISMATCH"
# res.message == "Experiment 'ENCSR000BST' targets TF 'GATA3', which does not match requested 'ESR1'."
```

---

## 7. Authentic ESR1 Metadata Audit & Status of Proposed Study

A live query against the official ENCODE portal identified the authentic released ESR1 ChIP-seq experiment in MCF-7 cells:
* **Authentic Experiment:** `ENCSR463GOT`
  - **Target:** `ESR1`
  - **Biosample:** `MCF-7`
  - **Assembly:** `GRCh38`
  - **Conservative IDR Peaks File:** `ENCFF185UUT` (Output type: `conservative IDR thresholded peaks`, MD5: `464964b5a099bdd901dfab66de064b21`)

### Definitive Operational Mandate:
1. **NOT AUTHORIZED FOR EXECUTION:** The *TFF1* / *ESR1* study remains **NOT AUTHORIZED FOR EXECUTION** during this release-preparation phase.
2. **Zero Ingestion / Execution:** No new datasets have been downloaded, and no biological experiments have been run or fabricated.
3. **Prospective Benchmark Firewall Protected:** The frozen prospective benchmark cohort locked in Phase 4C-1 (`ENCFF475LFH`, `ENCFF235PEX`, `ENCFF874PIO`, `ENCFF626VDA`, `ENCFF761OCJ`, `ENCFF572YTS`, `ENCFF148JKK`, `ENCFF635XCI`, `ENCFF171OJF`, `ENCFF771FRJ`) remains 100% firewalled, unblinded, and unexecuted.

---

## 8. Sign-off

**Lead Biological Data Auditor:** Computational Genomics Methodologist  
**Resolution Verdict:** **SECOND ENCODE EXPERIMENT ATTRIBUTION ERROR FORMALLY CORRECTED AND DOCUMENTED — ZERO BIOLOGICAL DATA INGESTION AUTHORIZED**
