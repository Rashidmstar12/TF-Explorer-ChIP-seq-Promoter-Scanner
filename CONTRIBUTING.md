# Contributing to TF-Explorer

Thank you for your interest in contributing to `TF-Explorer`! We welcome contributions from computational biologists, bioinformaticians, software engineers, and researchers.

## Code of Conduct

All contributors and participants are expected to adhere to the [Code of Conduct](CODE_OF_CONDUCT.md).

## How to Contribute

### Reporting Bugs and Provenance Discrepancies
If you encounter a bug, incorrect coordinate calculation, or dataset metadata misattribution:
1. Search existing issues to ensure the defect has not already been reported.
2. Open a new issue with a clear descriptive title.
3. Provide:
   - Python version and operating system.
   - Exact command line or Python code snippet.
   - Target gene, transcript ID, assembly (e.g. GRCh38), and cell line context.
   - Expected versus observed output.
   - Relevant error tracebacks or logs.

### Requesting Features
To suggest enhancements or new transcription factor binding data sources:
1. Open an issue describing the proposed feature.
2. Explain the scientific use case and biological rationale.
3. Suggest relevant public APIs or data formats if applicable.

### Developing and Submitting Pull Requests
1. **Fork and Clone:** Fork the repository on GitHub and clone your fork locally.
2. **Environment Setup:** Create a Python 3.10+ virtual environment and install dependencies:
   ```bash
   python -m venv venv
   source venv/bin/activate  # Or `venv\Scripts\activate` on Windows
   pip install -e ".[dev,test]"
   ```
3. **Coding Standards:**
   - Adhere to PEP 8 style conventions.
   - Write docstrings for all public classes, methods, and functions.
   - Ensure negative-strand promoter coordinates expand correctly into numerically higher coordinates.
   - Enforce dataset provenance using `tf_explorer.provenance`.
4. **Running Automated Tests:**
   All pull requests must pass the complete automated pytest suite:
   ```bash
   pytest -v
   ```
5. **Submitting PR:** Submit your pull request to the `main` branch with a concise explanation of your changes and test evidence.
