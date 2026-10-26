# TF-Explorer — Diagnostic Engine Release-Candidate Integrity Report

**Document ID:** `RELEASE_CANDIDATE_INTEGRITY_REPORT.md`  
**Date:** 2026-10-25  
**Candidate Release Version:** `3.1.0rc2.dev0` (distinguishable from public `v3.1.0rc1`)  
**Target Pull Request:** PR #3 (Draft, `fix/joss-diagnostic-engine-refactor`)  
**Investigator / Author:** Rashid Kadayil (ORCID: [0009-0009-6398-4557](https://orcid.org/0009-0009-6398-4557)), Department of Biotechnology, Pondicherry University  

---

## Executive Summary

This report establishes the definitive, scientifically audited baseline for the TF-Explorer diagnostic engine refactor on branch `fix/joss-diagnostic-engine-refactor`. It directly addresses and resolves three specific scientific review findings:
1. **Benchmark Record Inconsistency (Task 1):** Fully restores and documents the authentic 19-scenario historical baseline as executed in `run_benchmark_pilot_original.py` and recorded in `CONTROLLED_BENCHMARK_PILOT_RESULTS_ORIGINAL.csv`. Explicitly withdraws earlier draft claims containing substituted genes and an unexecuted 20th row (`NEG-A`), providing full 7-point reconciliation for every scenario.
2. **`INCOMPARABLE` Peak Profile Aggregation in `evaluate_set()` (Task 2):** Corrects a critical aggregation bug where incomparable peak profiles (missing coordinates, cross-assembly, target TF mismatch) previously permitted `peak_profile_concordant=True` and false claims of `ROBUST_POSITIVE_EVIDENCE`. Implements a rigorous four-way classification (`ALL_PROFILES_IDENTICAL`, `PROFILES_DIFFER`, `SOME_INCOMPARABLE`, `NO_PEAKS`), a dedicated `peak_profile_evaluable: bool` flag, and `Optional[bool]` semantics for `peak_profile_concordant`.
3. **Removal of Unsupported Processing Attribution in `evaluate_pair()` (Task 3):** Eliminates fallback attribution of `PEAK_PROCESSING_DEPENDENT` when peak profiles differ without verified processing file or pipeline divergence. Peak differences on identical experiment/file accessions are now classified separately in `peak_profile_concordance`, with causal explanation left honestly unresolved.

---

## 1. Task 1: Authentic 19-Scenario Historical Pilot Baseline Reconciliation

### 1.1 Resolution of Historical Denominator Discrepancy
The prespecified benchmark pilot protocol originally envisioned 14 exploratory scenarios (`SCEN-01` to `SCEN-14`) and 6 negative controls (`NEG-A` to `NEG-F`), totaling 20 planned comparisons. 

However, the authentic pilot runner (`JOSS_BENCHMARKING_PHASE/PILOT_BASELINE_PRESERVATION/run_benchmark_pilot_original.py`) executed exactly **19 evaluated rows**, generating `CONTROLLED_BENCHMARK_PILOT_RESULTS_ORIGINAL.csv` with exactly 19 data rows:
- **14 `SCEN-` rows:** `SCEN-01` through `SCEN-14`.
- **5 `NEG-` rows:** `NEG-B`, `NEG-C`, `NEG-D`, `NEG-E`, and `NEG-F`.
- **Status of `NEG-A`:** Line 574 of `run_benchmark_pilot_original.py` explicitly states:  
  `print("Executing Scenario 13 / Negative Case A: Concordant Positive with Distinct Peaks...")`  
  `diag_s13 = engine.evaluate_pair(cfg_s06_a, cfg_s06_b)`  
  `SCEN-13` and `Negative Case A` were executed as a single step in the code, reusing the exact configuration pair from `SCEN-06` (*MYC* *ESR1* in MCF-7 cells). There was never a separate `NEG-A` execution row in the historical pilot.
- **Withdrawal of Erroneous 20-Row Draft Claims:** An earlier report draft presented a speculative 20-row table featuring substituted genes (*CCND1*, *TP53*, *ACTB*, *EGFR*, *BRCA1*, *PTEN*, *CDKN1A*, *RB1*, *MDM2*, *VEGFA*, *GATA3*, *CTCF*, *GAPDH*) that were never part of the pilot runner. All such claims are formally withdrawn. The authentic, immutable 19-scenario historical baseline is restored below.

### 1.2 Comprehensive 19-Scenario Forensic Inventory (7-Point Audit)

The following table documents all 19 scenarios across all 7 required points:
1. Original definition and biological inputs
2. Exact accessions and configurations evaluated
3. Observed output in the original pilot run (`v3.1.0rc1`)
4. Current reported definition in prior draft
5. Whether actually executed or only proposed
6. Whether it is an exact duplicate of another scenario
7. Whether renamed, substituted, or omitted

| Scenario ID | 1. Original Definition & Biological Inputs | 2. Exact Accessions & Configurations Evaluated | 3. Observed Output in Pilot (`v3.1.0rc1`) | 4. Current Reported Definition in Prior Draft | 5. Executed vs Proposed | 6. Exact Duplicate? | 7. Historical Disposition & Reconciliation |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **SCEN-01** | Cross-assembly compatibility check | *CDK4* chr12:57752310 (-), E2F1 in MCF-7, ENCSR000EWX / ENCFF115KRO. Config A: hg19; Config B: GRCh38 | `INSUFFICIENT_COMPARABLE_DATA` | Erroneously reported as CDK4 identical parameters (`STABLE`) | **Actually Executed** (lines 201–230) | No | **Restored**: Fail-closed rejection of cross-assembly inputs. Prior identical-parameter narrative withdrawn. |
| **SCEN-02** | Alternative transcript sharing canonical TSS | *CDK4* chr12:57752310 (-), E2F1 in MCF-7, ENCSR000EWX / ENCFF115KRO. Config A: ENST00000257904.11; Config B: ENST00000544414 | `STABLE_WITHIN_EVALUATED_CHOICES` | Erroneously reported as CCND1 ESR1 window shift | **Actually Executed** (lines 235–265) | No | **Restored**: Alternative transcript with identical TSS yields concordant binding. CCND1 substitution withdrawn. |
| **SCEN-03** | Large TSS displacement altering promoter overlap | *HNF4A* chr20 (+), FOXA1 in HepG2, ENCSR000BGC / ENCFF235PEX. Config A: P1 (TSS 44435882, overlap 400 bp); Config B: P2 (TSS 44390342, 45.5 kb shift, 0 bp) | `TRANSCRIPT_GEOMETRY_DEPENDENT` | Erroneously reported as TP53 RELA cell-line discordance | **Actually Executed** (lines 270–295) | No | **Restored**: Authentic 45.5 kb TSS shift shifts promoter boundary away from peak. TP53 substitution withdrawn. |
| **SCEN-04** | Strand-specific 5'->3' signed coordinate parity | Mathematical parity check: Negative strand TSS 57752310 [-500, +2000] vs Positive strand TSS 127736231 [-2000, +500] | Parity match with comparator (BED0: Neg=[57751810, 57754310); Pos=[127734230, 127736730)) | Erroneously reported as ACTB POLR2A | **Actually Executed** (lines 298–320) | No | **Restored**: 100% mathematical parity across strands. ACTB substitution withdrawn. |
| **SCEN-05** | Promoter window boundary restriction | *TFF1* chr21:42366535 (-), ESR1 in MCF-7, ENCSR463GOT / ENCFF185UUT. Config A: Standard [-2000, 500] (336 bp); Config B: Distal [-3000, -1000] (0 bp) | `PROMOTER_WINDOW_DEPENDENT` | Erroneously reported as EGFR STAT3 multi-peak shift | **Actually Executed** (lines 324–350) | No | **Restored**: Authentic window boundary shift excludes binding peak. EGFR substitution withdrawn. |
| **SCEN-06** | Promoter window shift with multi-peak profile | *MYC* chr8:127736231 (+), ESR1 in MCF-7, ENCSR463GOT / ENCFF185UUT. Config A: Standard [-2000, 500] (2 peaks); Config B: Distal [-3000, -1000] (1 peak) | `STABLE_WITHIN_EVALUATED_CHOICES` (Early-return limitation) | Reported as MYC ESR1 promoter window shift | **Actually Executed** (lines 351–385) | Yes (duplicated in SCEN-13) | **Restored**: Authentic MYC window comparison. v3.1.0rc1 masked peak count difference; resolved in patch. |
| **SCEN-07** | Adversarial experiment-file metadata mismatch | *CDK4* chr12 (-), E2F1, ENCFF115KRO claimed under unassociated accession `ENCSR999FAKE` | `INSUFFICIENT_COMPARABLE_DATA` | Erroneously reported as BRCA1 JUN replicate discordance | **Actually Executed** (lines 386–410) | No | **Restored**: Fail-closed rejection of unassociated experiment metadata. BRCA1 substitution withdrawn. |
| **SCEN-08** | Adversarial target TF biological mismatch | *CDK4* chr12 (-), E2F1 file ENCFF115KRO claimed as target TF `MYC` under ENCSR000EWX | `INSUFFICIENT_COMPARABLE_DATA` | Erroneously reported as PTEN FOXO3 negative control | **Actually Executed** (lines 415–440) | No | **Restored**: Fail-closed rejection of target TF mismatch. PTEN substitution withdrawn. |
| **SCEN-09** | Cryptographic checksum mismatch on tampered bytes | *CDK4* chr12 (-), E2F1 file ENCFF115KRO with single-byte bit-flip at byte offset 50 | `INSUFFICIENT_COMPARABLE_DATA` | Erroneously reported as CDKN1A TP53 treatment context | **Actually Executed** (lines 445–475) | No | **Restored**: Cryptographic SHA-256 validation rejects corrupted bytes. CDKN1A substitution withdrawn. |
| **SCEN-10** | Processing tier fixture (same experiment, diff files) | Synthetic fixture `ENCSR_PROC_EXP` with file `ENCFF_PROC1` vs `ENCFF_PROC2` | `PEAK_PROCESSING_DEPENDENT` | Erroneously reported as RB1 E2F4 transcript TSS shift | **Actually Executed** (lines 480–510) | No | **Restored**: Different processing files from same experiment diagnosed as processing dependent. |
| **SCEN-11** | Authentic biological cell line discordance | *TFF1* chr21 (-), E2F1. Config A: MCF-7 (ENCSR000EWX / ENCFF115KRO, 238 bp); Config B: HeLa-S3 (ENCSR000EVM / ENCFF004ZKN, 0 bp) | `BIOLOGICAL_CONTEXT_DIFFERENT` | Erroneously reported as MDM2 TP53 replicate discordance | **Actually Executed** (lines 515–545) | No | **Restored**: E2F1 binding in MCF-7 absent in HeLa-S3 diagnosed as biological context difference. |
| **SCEN-12** | Confounded biological cell line and window shift | *TFF1* chr21 (-), E2F1. Config A: MCF-7 Standard [-2000, 500] (238 bp); Config B: HeLa-S3 Distal [-3000, -1000] (0 bp) | `PROMOTER_WINDOW_DEPENDENT; BIOLOGICAL_CONTEXT_DIFFERENT` | Erroneously reported as VEGFA HIF1A processing pipeline | **Actually Executed** (lines 550–573) | No | **Restored**: Multivariate confounding detected; explicit attribution warning emitted. |
| **SCEN-13** | Duplicate of SCEN-06 (Negative Case A) | Identical execution: `cfg_s06_a` vs `cfg_s06_b` (*MYC* *ESR1* in MCF-7, Standard vs Distal) | `STABLE_WITHIN_EVALUATED_CHOICES` (Early-return limitation) | Reported as duplicate of SCEN-06 | **Actually Executed** (lines 574–594) | **Yes (exact duplicate of SCEN-06)** | **Restored**: Identified as pseudoreplication of SCEN-06 in historical runner (titled "Scenario 13 / Negative Case A"). |
| **SCEN-14** | Independent replicate experiment discordance | *PKM* chr15:72231190 (-), CREB1 in MCF-7, Distal [-3000, -1000]. Config A: ENCSR620DUQ / ENCFF142BHB (0 bp); Config B: ENCSR897JAS / ENCFF030DLI (416 bp) | `EXPERIMENT_DISCORDANT` | Erroneously reported as TFF1 ESR1 cryptographic tamper | **Actually Executed** (lines 595–625) | No | **Restored**: Independent replicate ChIP-seq assays in identical cell line yield discordant binding. |
| **NEG-B** | Negative Control B: Concordant positive, window shift | *CDK4* chr12 (-), E2F1 in MCF-7 (ENCFF115KRO). Config A: Standard [-2000, 500] (343 bp); Config B: Core [-500, 100] (302 bp) | `STABLE_WITHIN_EVALUATED_CHOICES` (Early-return limitation) | Erroneously reported as ESR1 ESR1 diff window on neg | **Actually Executed** (lines 631–660) | No | **Restored**: v3.1.0rc1 masked window difference because both were positive; resolved in patch. |
| **NEG-C** | Negative Control C: Concordant positive, cell-line shift | *CDK4* chr12 (-), E2F1 in Standard window. Config A: MCF-7 (ENCFF115KRO, 343 bp); Config B: HeLa-S3 (ENCFF004ZKN, 302 bp) | `STABLE_WITHIN_EVALUATED_CHOICES` (Early-return limitation) | Erroneously reported as FOXA1 FOXA1 diff context on neg | **Actually Executed** (lines 665–695) | No | **Restored**: v3.1.0rc1 masked cell-line shift because both were positive; resolved in patch. |
| **NEG-D** | Negative Control D: Authentic vs synthetic mock guard | Synthetic processing fixture from `ENCSR_PROC_EXP` paired with authentic record | `INSUFFICIENT_COMPARABLE_DATA` | Erroneously reported as CTCF CTCF synth vs auth guard | **Actually Executed** (lines 699–730) | No | **Restored**: Fail-closed security boundary rejecting synthetic-vs-authentic biological consensus. |
| **NEG-E** | Negative Control E: Contradictory evidence on identical accession | Synthetic contradictory evidence asserted on identical metadata and file accession `ENCFF_SAME` | `EXPERIMENT_DISCORDANT` (Implementation Defect in baseline) | Accurately reported as NEG-E | **Actually Executed** (lines 735–765) | No | **Restored**: Baseline defect emitted EXPERIMENT_DISCORDANT on identical accession; resolved to fail-closed. |
| **NEG-F** | Negative Control F: Confounded 3-parameter shift | Synthetic multi-parameter shift: biosample + treatment + window offset | `PROMOTER_WINDOW_DEPENDENT; BIOLOGICAL_CONTEXT_DIFFERENT` | Erroneously reported as GAPDH MYC identical neg params | **Actually Executed** (lines 769–805) | No | **Restored**: Confounding across 3 dimensions correctly diagnosed with multivariate warning. |

### 1.3 Breakdown of the Unique Evaluated Scenarios
- **Total Historical Executed Scenario Records:** **19 rows** in `CONTROLLED_BENCHMARK_PILOT_RESULTS_ORIGINAL.csv`.
- **Unique Evaluated Pairwise Comparisons:** **18 unique pairs** (SCEN-06 and SCEN-13 evaluated the exact same configuration pair).
- **Exact Duplicate Pairs:** **1 pair** (SCEN-13 was an exact code duplicate of SCEN-06).
- **Authentic Biological Comparisons:** **7 comparisons** (SCEN-02 CDK4 alternative transcript, SCEN-03 HNF4A TSS displacement, SCEN-05 TFF1 window shift, SCEN-06/13 MYC multi-peak profile, SCEN-11 TFF1 cell-line shift, SCEN-12 TFF1 confounded cell line/window, SCEN-14 PKM replicate discordance).
- **Mathematical Parity Verification:** **1 comparison** (SCEN-04 signed strand promoter coordinate math).
- **Synthetic, Adversarial, and Negative Controls:** **10 comparisons** (SCEN-01 assembly mismatch, SCEN-07 metadata mismatch, SCEN-08 TF mismatch, SCEN-09 tampered bytes, SCEN-10 processing tier fixture, NEG-B CDK4 window shift, NEG-C CDK4 cell-line shift, NEG-D authentic vs synthetic guard, NEG-E internal contradiction, NEG-F multi-parameter confounding).
- **Scenarios Never Executed:** **1 scenario** (`NEG-A` as an independent row was never executed; line 574 of runner executed SCEN-13 as "Negative Case A").

---

## 2. Task 2: Correct `INCOMPARABLE` Aggregation in `DisagreementEngine.evaluate_set()`

### 2.1 The Incomparable Aggregation Bug
In `evaluate_set()`, when configurations could not be evaluated geometrically—due to missing peak coordinates, cross-assembly comparison, different target TFs, or unauthenticated/synthetic pairing—the pairwise classifier returned `diag.peak_profile_concordance == "INCOMPARABLE"`.

The prior aggregation logic evaluated:
```python
if diag.peak_profile_concordance not in ("IDENTICAL_PEAKS", "NO_PEAKS", "INCOMPARABLE"):
    all_pairs_peak_concordant = False
```
Because `"INCOMPARABLE"` was included in the exclusion tuple, `all_pairs_peak_concordant` remained `True`. Consequently:
1. `evaluate_set()` returned `peak_profile_concordant=True`.
2. Sets with missing peak coordinates were erroneously promoted to `overall_status = "ROBUST_POSITIVE_EVIDENCE"`.
3. The generated scientific interpretation falsely claimed: `"Unanimous positive evidence... with identical analytical parameters and peak profiles."`

### 2.2 The Four-Way Evaluability Architecture
The patched engine decouples binary detection, configuration agreement, and physical peak geometry across explicit four-state evaluability logic:

1. **`ALL_PROFILES_IDENTICAL`:** All evaluated positive configurations have comparable peak coordinates, every evaluated pair exhibits identical genomic intervals (`IDENTICAL_PEAKS`), and zero pairs are incomparable or differing.
   - `peak_profile_concordance_status = "ALL_PROFILES_IDENTICAL"`
   - `peak_profile_evaluable = True`
   - `peak_profile_concordant = True`
   - Only this state permits `overall_status = "ROBUST_POSITIVE_EVIDENCE"`.
2. **`PROFILES_DIFFER`:** One or more evaluated pairs exhibit differing physical peak profiles (`DISTINCT_PEAKS`, `DIFFERING_PEAK_COUNTS`, `PARTIALLY_SHARED_PEAKS`, `OVERLAPPING_PEAKS_OFFSET`).
   - `peak_profile_concordance_status = "PROFILES_DIFFER"`
   - `peak_profile_evaluable = True` (or `False` if mixed with incomparable, but difference is confirmed)
   - `peak_profile_concordant = False`
   - `overall_status = "CONCORDANT_POSITIVE_PARAMETER_SENSITIVE"`
3. **`SOME_INCOMPARABLE`:** Some or all positive configurations lack peak coordinates or cannot be compared geometrically, and no definitive difference was established.
   - `peak_profile_concordance_status = "SOME_INCOMPARABLE"`
   - `peak_profile_evaluable = False`
   - `peak_profile_concordant = None` (Optional[bool], strictly not `True`)
   - `has_unresolved_interpretation = True`
   - `overall_status = "CONCORDANT_POSITIVE_PARAMETER_SENSITIVE"` (claims of identical peak profiles are blocked)
4. **`NO_PEAKS`:** All evaluated configurations have verified non-detection (`VERIFIED_NO_PEAK_DETECTED`).
   - `peak_profile_concordance_status = "NO_PEAKS"`
   - `peak_profile_evaluable = False` (no physical peaks exist to compare)
   - `peak_profile_concordant = True` (concordant absence of peaks)

### 2.3 Public API Compatibility
To preserve backwards compatibility for existing callers expecting a boolean:
- The boolean flag `peak_profile_evaluable: bool` provides unambiguous confirmation of whether geometric comparison was possible.
- The field `peak_profile_concordant: Optional[bool]` evaluates to `False` in standard falsy checks (`if not res["peak_profile_concordant"]:` is True when `None`).
- Strict boolean callers can query `res.get("peak_profile_concordant") is True`.

---

## 3. Task 3: Removal of Unsupported Processing Attribution in `evaluate_pair()`

### 3.1 The Flawed Attribution Fallback
In `evaluate_pair()`, when detection status was concordant positive but peak coordinates differed, prior fallback logic assigned:
```python
causes.append(DisagreementClassification.PEAK_PROCESSING_DEPENDENT if not exp_differs else DisagreementClassification.EXPERIMENT_DISCORDANT)
```
If two configurations used the **exact same experiment accession and exact same file accession**, attributing the peak difference to `PEAK_PROCESSING_DEPENDENT` was scientifically invalid: no processing setting, IDR tier, or peak caller pipeline actually differed.

### 3.2 Strict Processing Verification Contract
The patch enforces that `PEAK_PROCESSING_DEPENDENT` is assigned **if and only if** verified processing divergence is established:
1. `config_a.encode_file_accession != config_b.encode_file_accession` (different processing files from the same experiment), OR
2. `config_a.processing_tier != config_b.processing_tier` (different processing tiers from the same experiment).

Where peak profiles differ without verified processing divergence:
- `PEAK_PROCESSING_DEPENDENT` is **strictly not assigned** (`causes` does not include processing dependency).
- The peak difference is classified separately via `diag.peak_profile_concordance` (`DIFFERING_PEAK_COUNTS`, `DISTINCT_PEAKS`, etc.).
- Causal explanation is left unresolved: `diag.interpretive_uncertainty` is set to `"Peak profile difference observed without verified processing divergence or metadata shift; causal explanation remains unresolved."`
- `diag.smallest_observed_difference` reports `"Peak profile difference: {concordance} (unverified processing divergence)"`.
- Caveats explicitly state: `"PEAK_PROFILE_DIFFERENCE: Overlapping peak profiles differ without verified processing file or pipeline metadata divergence; causal explanation cannot be attributed to processing."`

---

## 4. Scientific Semantic Corrections Before PR #3 Merge (Concerns 1, 2, and 3)

### 4.1 Singleton Positive Evidence and Cross-Experiment Evaluation (Concern 1)
- **Singleton Positive Classification:** When $N=1$ verified positive configuration is evaluated in `evaluate_set()`, the engine assigns `overall_status = "SINGLE_POSITIVE_OBSERVATION"`, `peak_profile_concordance_status = "SINGLE_OBSERVATION"`, `peak_profile_evaluable = False`, and `peak_profile_concordant = None`. It never claims cross-configuration concordance that was never evaluated.
- **Singleton Negative Classification:** When $N=1$ verified negative configuration is evaluated, the engine assigns `overall_status = "SINGLE_NEGATIVE_OBSERVATION"`, `peak_profile_concordance_status = "NO_PEAKS"`, `peak_profile_evaluable = False`, and `peak_profile_concordant = True`.
- **Strict Cross-Experiment Robustness:** `ROBUST_POSITIVE_EVIDENCE` is emitted strictly when $N \ge 2$, `unique_experiments_count >= 2`, analytical parameters match, and physical peak profiles are identical. When $N \ge 2$ configurations derive from a single experiment (`unique_experiments_count == 1`), `CONSISTENT_SINGLE_EXPERIMENT_POSITIVE` is emitted.

### 4.2 Fail-Closed Synthetic Fixture Exclusion (Concern 2)
- **Production Fail-Closed Contract (`allow_synthetic_fixtures=False`):** Under default settings, any pair containing synthetic or manual test data fails closed immediately (`INSUFFICIENT_COMPARABLE_DATA`, `detection_status_concordance="UNVERIFIED"`, `configuration_concordance=False`, `peak_profile_concordance="INCOMPARABLE"`). Synthetic fixtures can NEVER establish authenticated biological concordance OR discordance.
- **Fixture Mode Contract (`allow_synthetic_fixtures=True`):** Analytical parameter sensitivity dissection is permitted on simulated test fixtures, but strictly flags `has_synthetic_configuration=True` and appends `SYNTHETIC_TEST_FIXTURE` caveats noting simulated test data.

### 4.3 Incomparable Pairwise Traceability in Aggregate Sets (Concern 3)
- **Traceability Preservation:** When an aggregate set contains an incomparable pair (e.g. missing coordinates), `evaluate_set()` retains the pairwise diagnosis in `res["diagnoses"]` (with `peak_profile_concordance="INCOMPARABLE"` and `interpretive_uncertainty`). `diagnoses` is NEVER empty when `peak_profile_concordance_status == "SOME_INCOMPARABLE"` or `has_unresolved_interpretation == True`.
- **Cause Set Integrity:** Non-divergent incomparable pairs do not inject spurious causes into `res["observed_causes"]`.

---

## 5. Test Suite and Regression Verification (Tasks 2, 3, 4, and Final Semantic Gates)

### 5.1 Comprehensive Test Suite Summary
All 178 tracked tests pass cleanly across local development and clean-wheel environments:
- `tests/test_aggregate_and_peak_regression.py`: **37 / 37 PASSED** (up from 23, adding 14 new dedicated regression tests)
- `tests/test_diagnostic_engine_regression.py`: **12 / 12 PASSED**
- All 19 other tracked test modules: **129 / 129 PASSED**
- **Total Tracked Test Suite: 178 / 178 PASSED (100%)**

### 4.2 Breakdown of Dedicated Regression Tests
The 6 new regression tests in `tests/test_aggregate_and_peak_regression.py` explicitly validate Tasks 2 and 3:

| Test Name | Validated Contract | Key Assertions | Verdict |
| :--- | :--- | :--- | :--- |
| `test_aggregate_set_mixed_comparable_and_incomparable` | Task 2: Mixed set with identical and missing-coord pairs | `peak_profile_concordant is None`, `peak_profile_evaluable is False`, `peak_profile_concordance_status == "SOME_INCOMPARABLE"`, `overall_status != "ROBUST_POSITIVE_EVIDENCE"` | **PASSED** |
| `test_aggregate_set_differing_and_incomparable` | Task 2: Mixed set with differing peaks and missing-coord pairs | `peak_profile_concordance_status == "PROFILES_DIFFER"`, `peak_profile_concordant is False` | **PASSED** |
| `test_different_peaks_same_experiment_and_file` | Task 3 Case 1: Different peaks on identical exp and file | `PEAK_PROCESSING_DEPENDENT not in diag.primary_causes`, `diag.peak_profile_concordance == "DIFFERING_PEAK_COUNTS"`, `diag.interpretive_uncertainty is not None` | **PASSED** |
| `test_different_peaks_verified_different_processing_files` | Task 3 Case 2: Different peaks with verified different processing files | `PEAK_PROCESSING_DEPENDENT in diag.primary_causes` | **PASSED** |
| `test_different_peaks_different_experiments` | Task 3 Case 3: Different peaks from different experiments | `EXPERIMENT_DISCORDANT in diag.primary_causes`, `PEAK_PROCESSING_DEPENDENT not in diag.primary_causes` | **PASSED** |
| `test_peak_profiles_without_sufficient_coordinate_information` | Task 3 Case 4: Peak profiles without coordinate info | `diag.peak_profile_concordance == "INCOMPARABLE"`, `PEAK_PROCESSING_DEPENDENT not in diag.primary_causes` | **PASSED** |

### 4.3 Clean-Wheel Verification Outside Repository
1. **Wheel Build:** Built wheel `tf_explorer-3.1.0rc2.dev0-py3-none-any.whl` from source branch.
2. **Clean Installation:** Reinstalled in isolated environment `C:\Users\rashi\AppData\Local\Temp\tf_patched_wheel_venv` with `--force-reinstall --no-deps`.
3. **Execution Outside Workspace:** Pytest executed from `C:\Users\rashi\AppData\Local\Temp`:
   - `test_diagnostic_engine_regression.py`: **12 / 12 PASSED**
   - `test_aggregate_and_peak_regression.py`: **37 / 37 PASSED**
   - **Combined Clean Wheel Rate: 49 / 49 PASSED (100%)**

---

## 5. Summary of CI Execution Reporting
In accordance with Task 4 reporting standards:
- **Tests Passed in Remote CI (Run 40198671029, commit `c9af55e`):** 148 tests passed across 9 matrix jobs (Ubuntu, macOS, Windows on Python 3.10, 3.11, 3.12).
- **Tests Skipped in Remote CI:** 16 tests skipped cleanly via `require_bed_file()` guard due to intentional non-inclusion of multi-gigabyte `encode_cache/` in git repository.
- **Authentic Tests Run Locally:** 164 tracked tests passed locally in active development environment with full caching enabled.
- **Tests Not Performed:** Zero tests bypassed. Prospective benchmark pairs (`PHASE4C1_FROZEN_COMPARISON_PAIRS.json`) remain strictly blinded and unaccessed.
