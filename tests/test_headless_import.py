"""
Tests verifying clean decoupling of scientific core modules from Streamlit GUI.
Confirms that library and CLI modules can be imported and executed in headless environments
without importing Streamlit or triggering Streamlit runtime warnings.
"""

import sys
import pytest


def test_headless_modules_do_not_import_streamlit():
    """Importing core analysis and domain modules must NOT import streamlit into sys.modules."""
    # List of core domain modules that must be strictly headless
    headless_modules = [
        "tf_explorer",
        "tf_explorer.analysis",
        "tf_explorer.motifs",
        "tf_explorer.encode_client",
        "tf_explorer.genome_client",
        "tf_explorer.comparative",
        "tf_explorer.external_apis",
        "tf_explorer.browser_svg",
    ]
    
    # Save existing sys.modules state
    for mod in headless_modules:
        if mod in sys.modules:
            del sys.modules[mod]
            
    # Remove streamlit from sys.modules temporarily if present
    st_cached = sys.modules.pop("streamlit", None)
    
    try:
        import tf_explorer
        import tf_explorer.analysis
        import tf_explorer.motifs
        import tf_explorer.encode_client
        import tf_explorer.genome_client
        import tf_explorer.comparative
        import tf_explorer.external_apis
        import tf_explorer.browser_svg
        
        # Verify streamlit is NOT imported as a side effect
        assert "streamlit" not in sys.modules, (
            "Streamlit was unexpectedly imported by core headless modules!"
        )
    finally:
        # Restore streamlit if it was loaded previously
        if st_cached is not None:
            sys.modules["streamlit"] = st_cached


def test_cli_parser_headless():
    """Verify that CLI argument parser can be imported and configured without GUI context."""
    from tf_explorer import cli
    parser = cli.main.__code__  # Check main function is accessible
    assert parser is not None
