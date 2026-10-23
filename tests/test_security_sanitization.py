"""
Tests for input validation, sanitization, and security boundaries.
Verifies prevention of path traversal, shell injection metacharacters, and malformed inputs.
"""

import pytest
from tf_explorer.analysis import validate_gene_symbol, validate_genome_assembly


@pytest.mark.parametrize("valid_symbol", [
    "TP53", "EGFR", "MYC", "DNMT3B", "CDKN2A", "BRCA1", "miR-21"
])
def test_validate_gene_symbol_valid(valid_symbol):
    """Valid gene symbols must pass unchanged."""
    assert validate_gene_symbol(valid_symbol) == valid_symbol


@pytest.mark.parametrize("malformed_symbol", [
    "../../etc/passwd",
    "..\\windows\\system32",
    "TP53; rm -rf /",
    "EGFR && echo hacked",
    "MYC | cat",
    "BRCA1\0hidden",
    "CDKN2A\nBAD",
    "   ",
    "",
    12345,
    None
])
def test_validate_gene_symbol_rejections(malformed_symbol):
    """Path traversal and command metacharacters must raise ValueError."""
    with pytest.raises(ValueError):
        validate_gene_symbol(malformed_symbol)


@pytest.mark.parametrize("valid_input,expected", [
    ("hg38", "hg38"),
    ("GRCh38", "hg38"),
    ("HG38", "hg38"),
    ("hg19", "hg19"),
    ("GRCh37", "hg19"),
    ("HG19", "hg19"),
    ("", "hg38"),
    (None, "hg38")
])
def test_validate_genome_assembly_valid(valid_input, expected):
    """Supported genome assemblies must be recognized and normalized."""
    assert validate_genome_assembly(valid_input) == expected


@pytest.mark.parametrize("invalid_genome", [
    "mm10", "dm6", "sacCer3", "hg38; rm -rf", "invalid_assembly"
])
def test_validate_genome_assembly_invalid(invalid_genome):
    """Unsupported assemblies must raise ValueError."""
    with pytest.raises(ValueError):
        validate_genome_assembly(invalid_genome)
