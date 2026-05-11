# GraphBench Challenge — Plan d'implémentation complet

## Objectif

Créer un système complet de réfutation de conjectures en théorie des graphes, avec :
- **Partie 1** : heuristique de base (score v6 unique)
- **Partie 2** : heuristique FunSearch-inspired avec portefeuille de scores, mutations riches, et sélection automatique
- **Rapport PDF** : ~10 pages avec résultats, graphes, analyses comparatives
- **README + scripts** d'exécution

## Architecture générale

```
PRO OC /
├── graphbench_part1.py          # [EXISTING, to polish] Partie 1 - heuristique simple
├── graphbench_solver.py         # [EXISTING, to polish] Partie 2 - FunSearch-inspired  
├── benchmark.xlsx               # Données (fourni)
├── run_part1.bat                # Script d'exécution Part 1
├── run_part2.bat                # Script d'exécution Part 2
├── results/
│   ├── results_part1.csv        # Résultats Part 1
│   └── results_part2.csv        # Résultats Part 2
└── rapport/
    └── rapport_graphbench.md    # Rapport (10 pages)
```

## Proposed Changes

### Core Logic Files

#### [MODIFY] graphbench_part1.py
- Partie 1 : heuristique simple avec score v6 unique
- Mutations de base : add_edge, remove_edge, subdivide, add_vertex, etc.
- Export CSV complet
- Pas de FunSearch, pas de portfolio

#### [MODIFY] graphbench_solver.py
- Partie 2 : portefeuille de 9 fonctions de score (FunSearch-inspired)
- Mutations enrichies : densify_local, attach_star, attach_cycle, etc.
- Sélection automatique de la meilleure fonction de score par conjecture
- Support : connected, tree, claw_free, bipartite, planar
- Export CSV enrichi

### New Files

#### [NEW] run_part1.bat
Script Windows pour exécuter la Partie 1

#### [NEW] run_part2.bat
Script Windows pour exécuter la Partie 2

#### [NEW] rapport/rapport_graphbench.md
Rapport 10 pages contenant :
- Introduction, description du problème
- Architecture Part 1 (heuristique v6)
- Architecture Part 2 (FunSearch-inspired portfolio)
- Tableau des résultats par conjecture
- Analyse comparative Part1 vs Part2
- Graphes contre-exemples trouvés
- Lien GitHub
- Conclusion

## Verification Plan

### Automated Tests
- `python graphbench_part1.py --input benchmark.xlsx --output results/results_part1.csv --seconds 5 --max-order 14`
- `python graphbench_solver.py --input benchmark.xlsx --output results/results_part2.csv --seconds 5 --max-order 14`

### Manual Verification
- Vérifier que les CSV sont générés correctement
- Vérifier que les contre-exemples sont valides (graph6 valid, class valid, strict violation)
- Comparer les scores Part1 vs Part2
