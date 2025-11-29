@echo off
echo Running tf-explorer demo for PWWP2A...
echo.

REM Ensure we are using the correct python environment
py -m tf_explorer.cli --gene PWWP2A --tf-list "E2F1,YY1" --jaspar-ids "MA0024.3,MA0095.2" --bed-output --plot-track --out demo_results

echo.
echo Demo complete! Results are in the 'demo_results' directory.
pause
