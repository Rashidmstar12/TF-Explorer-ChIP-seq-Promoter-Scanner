# TF-Explorer — Diagnostic Engine Architectural Patch Summary

**Document ID:** `DIAGNOSTIC_ENGINE_PATCH_SUMMARY.md`  
**Date:** 2026-10-25  
**Candidate Release Version:** `3.1.0rc2.dev0` (Commit `fix/joss-diagnostic-engine-refactor`)  
**Target Pull Request:** PR #3 (Draft)  
**Author / Investigator:** Rashid Kadayil (ORCID: [0009-0009-6398-4557](https://orcid.org/0009-0009-6398-4557)), Department of Biotechnology, Pondicherry University  

---

## 1. Executive Summary and Rationale

The JOSS Controlled Benchmark Pilot evaluated published TF-Explorer `v3.1.0rc1` across the 19 prespecified historical scenarios in `CONTROLLED_BENCHMARK_PILOT_RESULTS_ORIGINAL.csv` against native bedtools and an Advanced Reference Workflow. The pilot revealed:
1. **13 Exact Agreements (68.4%):** Robust identification of authentic biological differences, strand-specific coordinate calculation, transcript TSS shifts, and cryptographic fail-closed rejection of tampered/unauthenticated genomic data.
2. **5 Design Limitations (26.3%):** A binary concordant-positive early-return mechanism (`STABLE_WITHIN_EVALUATED_CHOICES`) masked dimensional parameter differences (promoter window shifts, biological cell-line shifts, peak multiplicity shifts) whenever both configurations had detected peaks.
3. **1 Implementation Defect (5.3%):** In scenario `NEG-E`, when contradictory evidence statuses occurred with identical experiment accession (`ENCSR_SAME`), an empty candidate cause list caused unconditional fallback to `EXPERIMENT_DISCORDANT`, falsely alleging assay discordance across identical experiments.

Following independent scientific review of PR #3, three additional critical integrity corrections were implemented:
1. **Historical Pilot Baseline Reconciliation (Task 1):** Restoration of the authentic 19-scenario historical baseline from `CONTROLLED_BENCHMARK_PILOT_RESULTS_ORIGINAL.csv`, resolving the denominator discrepancy (pseudoreplication of *MYC* in SCEN-06 / SCEN-13, and clarification that NEG-A was never an independent executed row) and withdrawing speculative substituted gene narratives from earlier drafts.
2. **Incomparable Aggregation Fix (Task 2):** Prevention of false `peak_profile_concordant=True` on sets with incomparable peak profiles (missing coordinates, cross-assembly, target TF mismatch) in `evaluate_set()`. Introduction of four-way classification (`ALL_PROFILES_IDENTICAL`, `PROFILES_DIFFER`, `SOME_INCOMPARABLE`, `NO_PEAKS`), `peak_profile_evaluable: bool`, and `peak_profile_concordant: Optional[bool]`.
3. **Removal of Unsupported Processing Attribution (Task 3):** Strict enforcement that `PEAK_PROCESSING_DEPENDENT` is assigned only upon verified processing file or pipeline divergence. Peak differences on matching experiment/file accessions are classified separately in `peak_profile_concordance`, with causal explanation left honestly unresolved.

---

## 2. Taxonomy Extensions (`tf_explorer/disagreement_taxonomy.py`)

### 2.1 Standardized Classification Enum
In `DisagreementClassification`:
```python
class DisagreementClassification(str, Enum):
    TRANSCRIPT_GEOMETRY_DEPENDENT = "TRANSCRIPT_GEOMETRY_DEPENDENT"
    PROMOTER_WINDOW_DEPENDENT = "PROMOTER_WINDOW_DEPENDENT"
    EXPERIMENT_DISCORDANT = "EXPERIMENT_DISCORDANT"
    BIOLOGICAL_CONTEXT_DIFFERENT = "BIOLOGICAL_CONTEXT_DIFFERENT"
    PEAK_PROCESSING_DEPENDENT = "PEAK_PROCESSING_DEPENDENT"
    INSUFFICIENT_COMPARABLE_DATA = "INSUFFICIENT_COMPARABLE_DATA"
    STABLE_WITHIN_EVALUATED_CHOICES = "STABLE_WITHIN_EVALUATED_CHOICES"
    UNEXPLAINED_EVIDENCE_DISCORDANCE = "UNEXPLAINED_EVIDENCE_DISCORDANCE"
```

### 2.2 Decoupled Diagnostic Architecture
`DisagreementDiagnosis` provides dedicated fields decoupling binary detection, configuration agreement, and physical peak geometry:
- `detection_status_concordance: str` (`"BOTH_POSITIVE"`, `"BOTH_NEGATIVE"`, `"DISCORDANT"`, `"UNVERIFIED"`)
- `configuration_concordance: bool` (strict parameter agreement: TSS, window, cell line, treatment, genetic modification, processing tier)
- `peak_profile_concordance: str` (`"IDENTICAL_PEAKS"`, `"DISTINCT_PEAKS"`, `"DIFFERING_PEAK_COUNTS"`, `"PARTIALLY_SHARED_PEAKS"`, `"OVERLAPPING_PEAKS_OFFSET"`, `"NO_PEAKS"`, `"INCOMPARABLE"`)
- `primary_causes: List[DisagreementClassification]` (causes of divergence; contains `STABLE_WITHIN_EVALUATED_CHOICES` iff detection and parameters are concordant)
- `interpretive_uncertainty: Optional[str]` (explicitly records unresolved causal mechanisms)

---

## 3. Engine Redesign (`tf_explorer/disagreement_engine.py`)

### 3.1 Strict Causal Attribution for Peak Processing (Task 3)
Attribution to `PEAK_PROCESSING_DEPENDENT` is restricted strictly to verified file or tier divergence:
```python
processing_differs = same_exp and (file_differs or tier_differs)
if processing_differs:
    causes.append(DisagreementClassification.PEAK_PROCESSING_DEPENDENT)
```
Where peak profiles differ without verified processing divergence:
- `PEAK_PROCESSING_DEPENDENT` is **not** assigned.
- `diag.peak_profile_concordance` records the exact geometric shift (`DIFFERING_PEAK_COUNTS`, `DISTINCT_PEAKS`, etc.).
- `diag.interpretive_uncertainty` flags that causal explanation remains unresolved.

### 3.2 Incomparable Aggregation Fix in `evaluate_set()` (Task 2)
In `evaluate_set()`, peak profile concordance across configurations is classified into four mutually exclusive states:
1. `ALL_PROFILES_IDENTICAL`: All evaluated pairs with peaks have identical physical coordinates (`peak_profile_concordant=True`, `peak_profile_evaluable=True`).
2. `PROFILES_DIFFER`: Confirmed difference in physical peak profiles (`peak_profile_concordant=False`, `peak_profile_evaluable=True`).
3. `SOME_INCOMPARABLE`: Missing peak coordinates or metadata prevents unanimous geometric comparison (`peak_profile_concordant=None`, `peak_profile_evaluable=False`). Prevents false claims of `ROBUST_POSITIVE_EVIDENCE`.
4. `NO_PEAKS`: Unanimous absence of peaks across all configurations (`peak_profile_concordant=True`, `peak_profile_evaluable=False`).

---

## 4. Test Suite and Clean-Wheel Verification

1. **Tracked Core Test Suite:** **164 / 164 PASSED (100%)**
2. **Dedicated Regression Tests:**
   - `test_diagnostic_engine_regression.py`: **12 / 12 PASSED**
   - `test_aggregate_and_peak_regression.py`: **23 / 23 PASSED** (includes 6 new dedicated tests for Tasks 2 and 3)
3. **Clean-Wheel Verification Outside Repository:**
   - Wheel: `tf_explorer-3.1.0rc2.dev0-py3-none-any.whl`
   - Target Environment: `tf_patched_wheel_venv` (outside repository root)
   - Test Verdict: **35 / 35 PASSED (100%)**
4. **Remote CI Compatibility:**
   - 9 / 9 matrix jobs passed on GitHub Actions (Run `40198671029`)
   - Clean offline skipping via `require_bed_file()` preserves deterministic CI behavior without uncommitted multi-gigabyte cache dependencies.
