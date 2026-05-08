# GraphBench Challenge — Réfutation automatique de conjectures en théorie des graphes

## Prérequis

```bash
pip install networkx numpy pandas openpyxl
```

## Structure du projet

```
PRO OC hadi/
├── graphbench_part1.py     # Partie 1 : heuristique simple (score v6)
├── graphbench_solver.py    # Partie 2 : FunSearch-inspired (portefeuille de scores)
├── benchmark.xlsx          # Benchmark de conjectures
├── run_part1.bat           # Exécution Partie 1 (60s/conjecture)
├── run_part2.bat           # Exécution Partie 2 (60s/conjecture)
├── results/                # Résultats CSV générés
└── rapport/                # Rapport PDF (~10 pages)
```

## Exécution rapide

### Partie 1 (heuristique de base)
```bash
python graphbench_part1.py --input benchmark.xlsx --output results/results_part1.csv --seconds 60 --max-order 18
```

### Partie 2 (FunSearch-inspired)
```bash
python graphbench_solver.py --input benchmark.xlsx --output results/results_part2.csv --seconds 60 --max-order 18
```

### Test rapide (5 secondes par conjecture)
```bash
python graphbench_part1.py --input benchmark.xlsx --output results/test_part1.csv --seconds 5 --max-order 14
python graphbench_solver.py --input benchmark.xlsx --output results/test_part2.csv --seconds 5 --max-order 14
```

## Description

### Partie 1 — Heuristique simple

- **Banque initiale** : graphes atlas (≤7 sommets), chemins, cycles, cliques, étoiles, roues, bipartis complets, arbres de Prüfer, graphes aléatoires Gnp
- **Mutations** : `add_edge`, `remove_edge`, `toggle_edge`, `add_vertex`, `remove_vertex`, `subdivide`, `add_path`, `add_clique`, `add_twins`, `line_graph`
- **Réparation** : connexité, arbre, claw-free
- **Score** : `score_v6_final_best` — marge normalisée + direction + structure + signal (y−f(x)) − complexité

### Partie 2 — FunSearch-inspired

- **Portefeuille de 9 fonctions de score** évaluées à chaque étape :
  - `v0` violation pure, `v1` marge normalisée, `v2` structure guidée
  - `v3` direction de conjecture, `v4` spécialisée par classe
  - `v5` invariants difficiles guidés, `v6` sparse/long, `v7` dense/clique, `v8` mixte final
- **Sélection automatique** : la meilleure valeur parmi les 9 est retenue à chaque évaluation
- **Mutations enrichies** : + `densify_local`, `sparsify_local`, `attach_star`, `attach_cycle`, `double_subdivide`
- **Classes supportées** : `connected`, `tree`, `claw_free`, `bipartite`, `planar`
- **Population évolutive** : top-140 candidats maintenus, mutations multi-étapes (1–3 steps)

## Invariants calculés

| Alias | Invariant |
|-------|-----------|
| `n` | order (nombre de sommets) |
| `m` | size (nombre d'arêtes) |
| `diam` | diameter |
| `rad` | radius |
| `delta` | minimum_degree |
| `Delta` | maximum_degree |
| `avg` | average_degree |
| `t` | triangle_number |
| `omega` | clique_number |
| `gamma` | domination_number |
| `gamma_t` | total_domination_number |
| `alpha` | independence_number |
| `tau` | vertex_cover_number |
| `mu` | matching_number |
| `kappa` | vertex_connectivity |
| `kappa_prime` | edge_connectivity |

Plus : `randic_index`, `harmonic_index`, `first_zagreb_index`, `second_zagreb_index`, `proximity`, `remoteness`, `largest_eigenvalue`, `largest_distance_eigenvalue`, `second_smallest_laplace_eigenvalue`

## Format de sortie CSV

| Colonne | Description |
|---------|-------------|
| `conjecture_id` | ID de la conjecture |
| `status` | FOUND / NOT_FOUND |
| `validation_status` | VALID_COUNTEREXAMPLE / ... |
| `score_variant` | Nom de la fonction de score utilisée |
| `counterexample_g6` | Graphe en format graph6 |
| `order`, `size` | Taille du contre-exemple |
| `counter_x`, `counter_y` | Valeurs des invariants |
| `violation_margin` | Marge de violation (> 0 si contre-exemple) |
| `cost` | Temps trouvé (ou 120 si non trouvé) |
| `execution_time_seconds` | Temps d'exécution réel |

## Scoring

- Si contre-exemple trouvé en `t` secondes → coût = `t`
- Si non trouvé en 60 secondes → pénalité = `120` secondes
- **Score total** = somme des coûts sur toutes les conjectures (plus bas = meilleur)

## Notes importantes

- Les contre-exemples fournis dans `benchmark.xlsx` ne sont **pas utilisés** par l'algorithme
- Le code fonctionne avec n'importe quel fichier benchmark de même format
- Les invariants exacts sont calculés pour n ≤ 20, approchés par greedy sinon
