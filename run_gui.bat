@echo off
echo Starting tf-explorer GUI...
echo.

REM Ensure we are using the correct python environment
py -m streamlit run tf_explorer/app.py

pause
