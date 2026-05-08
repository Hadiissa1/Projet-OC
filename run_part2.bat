@echo off
echo ============================================
echo  GraphBench - Partie 2 (FunSearch-inspired)
echo ============================================
mkdir results 2>nul

REM Detecter le Python disponible (venv en priorite)
set PYTHON=python
if exist "%~dp0.venv\Scripts\python.exe" set PYTHON=%~dp0.venv\Scripts\python.exe

echo Python utilise : %PYTHON%
echo.

%PYTHON% "%~dp0graphbench_solver.py" --input "%~dp0benchmark.xlsx" --output "%~dp0results\results_part2.csv" --seconds 60 --max-order 18

echo.
echo Termine ! Resultats dans results\results_part2.csv
pause
