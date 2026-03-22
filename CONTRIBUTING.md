# Contributing to TF-Explorer

Thank you for considering contributing to TF-Explorer! 🎉

All contributions are welcome: bug reports, feature requests, documentation improvements, and code patches.

## How to Contribute

### 1 · Report a Bug

Open an issue on [GitHub Issues](https://github.com/Rashidmstar12/TF-Explorer-ChIP-seq-Promoter-Scanner/issues) and include:

- A clear description of the problem.
- Steps to reproduce it (gene name, TFs used, command run).
- Your OS, Python version, and the version of TF-Explorer.
- Relevant error messages or screenshots.

### 2 · Suggest a Feature

Open an issue with the label **enhancement** and describe:

- What you want the tool to do.
- Why it would be useful.
- (Optional) A rough idea of how it could be implemented.

### 3 · Submit a Pull Request

1. **Fork** the repository and create a feature branch:

   ```bash
   git checkout -b feature/my-new-feature
   ```

2. **Install** the development dependencies:

   ```bash
   pip install -r requirements.txt
   pip install -e .
   ```

3. **Make your changes**, keeping them focused and well-documented.

4. **Run the tests** to make sure nothing is broken:

   ```bash
   python -m pytest test_simple.py test_metrics_logic.py test_comparative.py test_cell_line_comparison.py -v
   ```

5. **Commit** with a clear message and **push** your branch, then open a pull request against `main`.

## Code Style

- Follow [PEP 8](https://pep8.org/) for Python code.
- Use descriptive variable names.
- Add docstrings to new functions and classes.
- Keep functions focused and reasonably short.

## Questions?

Feel free to open an issue or start a discussion on GitHub.
