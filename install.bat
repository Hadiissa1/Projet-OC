@echo off
REM Script d'installation du projet GraphBench
echo Installation du projet GraphBench...

REM Copier les scripts Python depuis le dossier "projet oc"
set SRC=C:\Users\Hadi\OneDrive\Desktop\M1 Miage\M1 S2\OC\projet oc
set DST=C:\Users\Hadi\OneDrive\Desktop\M1 Miage\M1 S2\OC\PRO OC hadi

copy "%SRC%\graphbench_part1.py" "%DST%\graphbench_part1.py"
copy "%SRC%\graphbench_solver.py" "%DST%\graphbench_solver.py"
copy "%SRC%\benchmark.xlsx" "%DST%\benchmark.xlsx" 2>nul

REM Créer le dossier results
mkdir "%DST%\results" 2>nul

echo.
echo Verification des dependances Python...
python -c "import networkx, numpy, pandas, openpyxl; print('Toutes les dependances sont installees.')"
if errorlevel 1 (
    echo Installation des dependances...
    pip install networkx numpy pandas openpyxl
)

echo.
echo Installation terminee !
echo.
echo Pour lancer la Partie 1 :
echo   python "%DST%\graphbench_part1.py" --input "%DST%\benchmark.xlsx" --output "%DST%\results\results_part1.csv" --seconds 60 --max-order 18
echo.
echo Pour lancer la Partie 2 :
echo   python "%DST%\graphbench_solver.py" --input "%DST%\benchmark.xlsx" --output "%DST%\results\results_part2.csv" --seconds 60 --max-order 18
echo.
pause
