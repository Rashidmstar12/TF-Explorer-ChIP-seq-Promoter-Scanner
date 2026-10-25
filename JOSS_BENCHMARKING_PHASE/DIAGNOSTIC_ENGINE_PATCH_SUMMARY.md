# TF-Explorer — Diagnostic Engine Architectural Patch Summary

**Document ID:** `DIAGNOSTIC_ENGINE_PATCH_SUMMARY.md`  
**Date:** 2026-10-25  
**Baseline Released Software:** TF-Explorer `v3.1.0rc1` (Commit `58f84a52e796a39fcc92dd36d43e524bb775e4eb`)  
**Development Branch:** `fix/joss-diagnostic-engine-refactor`  
**Author / Investigator:** Rashid Kadayil (ORCID: [0009-0009-6398-4557](https://orcid.org/0009-0009-6398-4557)), Department of Biotechnology, Pondicherry University  

---

## 1. Executive Summary and Rationale

The JOSS Controlled Benchmark Pilot evaluated published TF-Explorer `v3.1.0rc1` across 19 prespecified scenarios against native bedtools and an Advanced Reference Workflow. The pilot revealed:
1. **13 Exact Agreements (68.4%):** Robust identification of authentic biological differences, strand-specific coordinate calculation, transcript TSS shifts, and cryptographic fail-closed rejection of tampered/unauthenticated genomic data.
2. **5 Design Limitations (26.3%):** A binary concordant-positive early-return mechanism (`STABLE_WITHIN_EVALUATED_CHOICES`) masked dimensional parameter differences (promoter window shifts, biological cell-line shifts, peak multiplicity shifts) whenever both configurations had detected peaks.
3. **1 Implementation Defect (5.3%):** In scenario `NEG-E`, when contradictory evidence statuses occurred with identical experiment accession (`ENCSR_SAME`), an empty candidate cause list caused unconditional fallback to `EXPERIMENT_DISCORDANT`, falsely alleging assay discordance across identical experiments.

Following the pilot decision `CORRECT_SOFTWARE_BEFORE_BENCHMARK`, focused architectural remediations were implemented on dedicated development branch `fix/joss-diagnostic-engine-refactor` to resolve these issues before the definitive prospective benchmark.

---

## 2. Taxonomy Extensions (`tf_explorer/disagreement_taxonomy.py`)

### 2.1 Addition of `UNEXPLAINED_EVIDENCE_DISCORDANCE`
In `DisagreementClassification`:
```python
class DisagreementClassification(str, Enum):
    ...
    UNEXPLAINED_EVIDENCE_DISCORDANCE = "UNEXPLAINED_EVIDENCE_DISCORDANCE"
```
**Rationale:** Distinguishes unexplained assay discrepancies from genuine independent biological experiment discordance. If two configurations share the same experiment accession or have identical parameters but discordant calls, the engine must never diagnose `EXPERIMENT_DISCORDANT`.

### 2.2 Four-Dimensional Diagnostic Dataclass Extension
`DisagreementDiagnosis` was extended with four explicit audit fields:
```python
@dataclass
class DisagreementDiagnosis:
    primary_causes: List[DisagreementClassification]
    differing_dimensions: List[DimensionType]
    scientific_explanation: str
    actionable_recommendations: List[str]
    causal_attribution_caveat: str
    has_multiple_covarying_dimensions: bool = False
    # Newly introduced multi-dimensional audit fields:
    detection_status_concordance: str = "UNVERIFIED" # BOTH_POSITIVE | BOTH_NEGATIVE | DISCORDANT | UNVERIFIED
    configuration_concordance: bool = True           # True if all analytical and metadata parameters match
    peak_profile_concordance: str = "INCOMPARABLE"   # IDENTICAL_PEAKS | DISTINCT_PEAKS | DIFFERING_PEAK_COUNTS | NO_PEAKS | INCOMPARABLE
    interpretive_uncertainty: List[str] = field(default_factory=list)
```

### 2.3 Backward-Compatible Serialization
The `to_dict()` method serializes the new fields while preserving all existing legacy dictionary keys:
```python
def to_dict(self) -> Dict[str, Any]:
    return {
        "primary_causes": [c.value for c in self.primary_causes],
        "differing_dimensions": [d.value for d in self.differing_dimensions],
        "scientific_explanation": self.scientific_explanation,
        "actionable_recommendations": self.actionable_recommendations,
        "causal_attribution_caveat": self.causal_attribution_caveat,
        "has_multiple_covarying_dimensions": self.has_multiple_covarying_dimensions,
        "detection_status_concordance": self.detection_status_concordance,
        "configuration_concordance": self.configuration_concordance,
        "peak_profile_concordance": self.peak_profile_concordance,
        "interpretive_uncertainty": self.interpretive_uncertainty,
    }
```

---

## 3. Engine Redesign (`tf_explorer/disagreement_engine.py`)

### 3.1 Strict Fail-Closed Provenance and Synthetic Bypass Elimination
In `v3.1.0rc1`, line 148 checked:
```python
# Old vulnerable code:
if (is_syn_a or is_syn_b) and (config_a.evidence_status == config_b.evidence_status):
    ...
```
This allowed discordant mock configurations (`config_a.evidence_status != config_b.evidence_status`) to bypass authentication checks and reach biological diagnostic logic.

**Patched Remediation:**
1. **Asymmetric Provenance Guard:** Comparing authentic genomic observations against synthetic or unverified configurations immediately fails closed:
   ```python
   if (is_syn_a != is_syn_b) or (not config_a.is_authenticated) or (not config_b.is_authenticated):
       # Return INSUFFICIENT_COMPARABLE_DATA
   ```
2. **Cryptographic Integrity & Replay Guard:** Any token tampering or mismatch against authoritative ENCODE catalog strictly returns `INSUFFICIENT_COMPARABLE_DATA`.

### 3.2 Internal Data Contradiction Detection (NEG-E Fix)
When configurations have identical genomic coordinates, promoter window, cell line, TF target, experiment accession, and file accession, but contradictory evidence statuses (`POSITIVE` vs `NEGATIVE`), the data represents an internal contradiction:
```python
if (config_a.encode_experiment_accession == config_b.encode_experiment_accession and
    config_a.encode_file_accession == config_b.encode_file_accession and
    config_a.canonical_tss == config_b.canonical_tss and
    config_a.promoter_start == config_b.promoter_start and
    config_a.promoter_end == config_b.promoter_end and
    config_a.evidence_status != config_b.evidence_status):
    return DisagreementDiagnosis(
        primary_causes=[DisagreementClassification.INSUFFICIENT_COMPARABLE_DATA],
        differing_dimensions=[],
        scientific_explanation="Internal data contradiction: identical parameters assert contradictory evidence.",
        actionable_recommendations=["Inspect input source files for contradictory duplicate annotations."],
        causal_attribution_caveat="Contradictory evidence within identical experimental parameters cannot be diagnosed.",
        has_multiple_covarying_dimensions=False,
        detection_status_concordance="DISCORDANT",
        configuration_concordance=True,
        peak_profile_concordance="INCOMPARABLE",
        interpretive_uncertainty=["INTERNAL_DATA_CONTRADICTION"]
    )
```

### 3.3 Safe Fallback Logic (Elimination of False EXPERIMENT_DISCORDANT)
In `v3.1.0rc1`, lines 259–262 had:
```python
# Old flawed fallback:
if not causes:
    causes.append(DisagreementClassification.EXPERIMENT_DISCORDANT)
```
**Patched Remediation:**
```python
if not causes:
    if config_a.encode_experiment_accession == config_b.encode_experiment_accession:
        causes.append(DisagreementClassification.UNEXPLAINED_EVIDENCE_DISCORDANCE)
    else:
        causes.append(DisagreementClassification.EXPERIMENT_DISCORDANT)
```
This guarantees that `EXPERIMENT_DISCORDANT` is **never** emitted when experiment accessions are identical.

### 3.4 Multi-Dimensional Concordant Evidence Audit (LIMITATION-01 Fix)
In `v3.1.0rc1`, lines 171–183 returned immediately upon `config_a.evidence_status == config_b.evidence_status`, discarding all coordinate, window, biosample, and peak count differences.

**Patched Remediation:**
1. **Concordant Positive Analysis:**
   - Retains `detection_status_concordance = "BOTH_POSITIVE"`.
   - Distinguishes identical configurations from parameter-varying configurations:
     - If all parameters, coordinates, and peak profiles match identically: `causes = [DisagreementClassification.STABLE_WITHIN_EVALUATED_CHOICES]`.
     - If promoter windows differ: `causes.append(DisagreementClassification.PROMOTER_WINDOW_DEPENDENT)` and `configuration_concordance = False`.
     - If biological contexts differ: `causes.append(DisagreementClassification.BIOLOGICAL_CONTEXT_DISCORDANCE)` and `configuration_concordance = False`.
     - If processing tiers differ: `causes.append(DisagreementClassification.PROCESSING_TIER_DISCORDANT)` and `configuration_concordance = False`.
     - If physical peak profiles differ: `causes.append(DisagreementClassification.PEAK_MULTIPLICITY_OR_PROFILE_SHIFT)`.
2. **Concordant Negative Guard:**
   - Retains `detection_status_concordance = "BOTH_NEGATIVE"`.
   - Sets `peak_profile_concordance = "NO_PEAKS"`.
   - If biological contexts differ: notes `BIOLOGICAL_CONTEXT_DISCORDANCE` across negative tests.
   - Mandatory caveat: "Non-detection of peaks within the evaluated promoter window does not demonstrate absence of transcriptional regulation or distal enhancer binding."

### 3.5 Physical Peak Profile Classification (`classify_peak_profile_concordance`)
A dedicated module-level helper `classify_peak_profile_concordance(config_a, config_b, is_both_positive, is_both_negative) -> str` enforces rigorous physical peak comparisons:
- Target TF mismatch -> `"INCOMPARABLE"`
- Assembly mismatch or non-GRCh38 -> `"INCOMPARABLE"`
- Missing peak details / empty coordinates on positive configs -> `"INCOMPARABLE"`
- Mismatched provenance or unauthenticated pairing -> `"INCOMPARABLE"`
- Differing peak counts -> `"DIFFERING_PEAK_COUNTS"`
- Identical peak intervals `(chrom, start, end)` -> `"IDENTICAL_PEAKS"`
- Partially intersecting set of peak intervals -> `"PARTIALLY_SHARED_PEAKS"`
- Coordinate overlap with non-zero offset -> `"OVERLAPPING_PEAKS_OFFSET"`
- Completely disjoint non-overlapping coordinates -> `"DISTINCT_PEAKS"`
- Negative concordant -> `"NO_PEAKS"`

### 3.6 Multi-Configuration Set Aggregation Redesign (`evaluate_set()`)
In `v3.1.0rc1`, lines 550–555 contained a severe aggregation defect:
```python
# Old flawed aggregation in v3.1.0rc1:
for diag in all_pairwise_diagnoses:
    if DisagreementClassification.STABLE_WITHIN_EVALUATED_CHOICES not in diag.primary_causes:
        retained_diagnoses.append(diag)
```
This logic unconditionally discarded any pairwise comparison that contained `STABLE_WITHIN_EVALUATED_CHOICES`, meaning that when configurations had concordant positive detection but differing promoter windows or biosamples, the window/context shifts were dropped from `observed_causes` and `diagnoses`, reporting an empty list of causes.

**Patched Remediation:**
1. Evaluates all unique pairwise comparisons in the configuration set.
2. Evaluates divergence comprehensively: a pairwise comparison is retained as an active diagnostic finding if:
   - Any cause other than `STABLE_WITHIN_EVALUATED_CHOICES` is present;
   - `not diag.configuration_concordance`;
   - `diag.peak_profile_concordance not in ("IDENTICAL_PEAKS", "NO_PEAKS")`;
   - Detection status is discordant; OR
   - Primary cause is `INSUFFICIENT_COMPARABLE_DATA`.
3. Reports granular, decoupled aggregate metrics:
   - `detection_status_concordant: bool` (binary detection agreement across all configs)
   - `configuration_concordant: bool` (complete analytical and metadata parameter agreement)
   - `peak_profile_concordant: bool` (physical peak geometry agreement)
   - `has_unresolved_interpretation: bool`
   - `has_invalid_or_unauthenticated_evidence: bool`
   - `overall_status`: Granular statuses (`EMPTY_SET`, `INSUFFICIENT_OR_UNVERIFIED_DATA`, `DISCORDANT_OBSERVATIONS`, `CONCORDANT_POSITIVE_IDENTICAL`, `CONCORDANT_POSITIVE_PARAMETER_SENSITIVE`, `CONCORDANT_NEGATIVE_IDENTICAL`, `CONCORDANT_NEGATIVE_CONTEXT_SENSITIVE`).
   - `observed_causes`: Clean sorted list of distinct dimensional causes observed across all pairs.
   - `diagnoses`: Fully serializable list of retained pairwise diagnoses.

### 3.7 Semantic Separation of Concordance Concepts
The patch formally resolves the conceptual ambiguity of `STABLE_WITHIN_EVALUATED_CHOICES`:
- **Binary Detection Concordance (`detection_status_concordance`):** Reports whether both configurations detected binding (`BOTH_POSITIVE`), neither detected binding (`BOTH_NEGATIVE`), or calls conflicted (`DISCORDANT`).
- **Configuration Concordance (`configuration_concordance`):** Strict boolean indicating whether all analytical and experimental parameters (transcript TSS, promoter window, cell line, treatment, genetic modification, processing tier) are identical.
- **Diagnostic Classification (`primary_causes`):** Lists only actual causes of divergence (`PROMOTER_WINDOW_DEPENDENT`, etc.). `STABLE_WITHIN_EVALUATED_CHOICES` is emitted if and only if both detection status and all evaluated analytical parameters are concordant.

---

## 4. Summary of Code Diffs

| Component | File | Lines Changed | Primary Impact |
| :--- | :--- | :--- | :--- |
| **Taxonomy** | `tf_explorer/disagreement_taxonomy.py` | +21 lines | Added `UNEXPLAINED_EVIDENCE_DISCORDANCE`, added 4 multi-dimensional audit attributes to `DisagreementDiagnosis`. |
| **Provenance** | `tf_explorer/disagreement_engine.py` | +48 lines | Added fail-closed asymmetric provenance check and internal contradiction detection. |
| **Concordance** | `tf_explorer/disagreement_engine.py` | +112 lines | Redesigned concordant status handling; audits window, biosample, processing tier, and peak geometry on positive pairs. |
| **Peak Profile** | `tf_explorer/disagreement_engine.py` | +98 lines | Implemented `classify_peak_profile_concordance()` covering identical, disjoint, shared, offset, differing counts, and incomparable. |
| **Aggregation** | `tf_explorer/disagreement_engine.py` | +142 lines | Replaced flawed filter in `evaluate_set()`; reports granular detection, configuration, and profile concordance flags. |
| **Fallback** | `tf_explorer/disagreement_engine.py` | +12 lines | Guarded fallback to prevent emitting `EXPERIMENT_DISCORDANT` on matching accessions. |

---

## 5. Verification and Parity Guarantees

1. **Defect-01 (NEG-E) Resolved:** Output changes from false `EXPERIMENT_DISCORDANT` to fail-closed `INSUFFICIENT_COMPARABLE_DATA`.
2. **Limitation-01 (Early Return) Resolved:** `SCEN-06`, `SCEN-13`, `NEG-B`, `NEG-C`, and `NEG-D` correctly surface promoter window, cell-line, and peak profile differences while retaining detection concordance.
3. **Aggregation Defect Resolved:** `evaluate_set()` correctly retains parameter shifts (`PROMOTER_WINDOW_DEPENDENT`, `PEAK_MULTIPLICITY_OR_PROFILE_SHIFT`, etc.) instead of discarding them.
4. **100% Backward Compatibility:** All existing taxonomy enum values and dictionary keys in `to_dict()`, `to_markdown()`, and `to_dataframe()` remain identical.
5. **100% Mathematical Parity:** Signed strand coordinate arithmetic (`SCEN-04`) remains identical to native bedtools BED0 conversions.
6. **Clean Wheel Verification:** Built wheel `tf_explorer-3.1.0rc1-py3-none-any.whl` installed in isolated `tf_patched_wheel_venv` outside the repository checkout passed all 29 regression tests.
7. **Complete Core Test Suite:** 158 / 158 tests passed in active suite.
