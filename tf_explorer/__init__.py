"""
TF-Explorer: Epigenetic Promoter Scanner & Regulatory Disagreement Engine.
"""

__version__ = "3.1.0rc1"

# Core public API exports
from tf_explorer.analysis import calculate_promoter_coordinates
from tf_explorer.candidate_screener import CandidateScreener
from tf_explorer.disagreement_engine import DisagreementEngine
from tf_explorer.disagreement_taxonomy import (
    DisagreementClassification,
    DisagreementDiagnosis,
    EvidenceConfiguration,
    EvidenceInterpretation,
    DimensionType,
)

__all__ = [
    "__version__",
    "calculate_promoter_coordinates",
    "CandidateScreener",
    "DisagreementEngine",
    "DisagreementClassification",
    "DisagreementDiagnosis",
    "EvidenceConfiguration",
    "EvidenceInterpretation",
    "DimensionType",
]
