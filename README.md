# Projet GraphBench Challenge

## Réfutation automatique de conjectures en théorie des graphes

Ce projet a été réalisé dans le cadre du module **Optimisation Combinatoire — Master 1 MIAGE**.

L’objectif du projet est de développer un programme capable de rechercher automatiquement des contre-exemples à des conjectures en théorie des graphes. Pour chaque conjecture, le programme génère des graphes candidats, calcule leurs invariants, vérifie si l’inégalité est strictement violée et exporte les résultats dans un fichier CSV.

---

## Contenu du dossier

Le dossier du projet est organisé de la manière suivante :

```text
PROJET-OC/
├── results/
│   ├── test_part1.csv
│   └── test_part2.csv
├── benchmark.xlsx
├── graphbench_part1.py
├── graphbench_solver.py
├── README.md
└── requirements.txt
```

### Description des fichiers

- `benchmark.xlsx` : fichier contenant les 100 conjectures à tester.
- `graphbench_part1.py` : script Python de la Partie 1, basé sur une heuristique simple de recherche locale.
- `graphbench_solver.py` : script Python de la Partie 2, version améliorée inspirée de FunSearch.
- `results/test_part1.csv` : fichier CSV contenant les résultats de la Partie 1.
- `results/test_part2.csv` : fichier CSV contenant les résultats de la Partie 2.
- `README.md` : documentation du projet.
- `requirements.txt` : liste des dépendances Python nécessaires.

---

## Principe général

Le programme lit le fichier `benchmark.xlsx`, reconstruit chaque conjecture, génère des graphes candidats, calcule les invariants de graphes nécessaires et vérifie automatiquement si l’inégalité associée à la conjecture est violée.

Un contre-exemple est accepté uniquement si :

- le graphe appartient à la classe demandée ;
- les invariants nécessaires sont calculés ;
- l’inégalité est strictement violée ;
- le graphe est exporté au format `graph6`.

Point important : la colonne `Counter example (g6)` éventuellement présente dans le benchmark n’est pas utilisée comme réponse. Les contre-exemples présents dans les fichiers CSV sont générés automatiquement par les scripts.

---

## Partie 1 — Heuristique simple

Le fichier `graphbench_part1.py` correspond à la première approche.

Cette version utilise :

- une banque initiale de graphes ;
- une fonction de score ;
- des mutations locales ;
- un mécanisme de réparation ;
- une population de graphes candidats ;
- une validation automatique des contre-exemples.

Commande d’exécution :

```bash
python graphbench_part1.py --input benchmark.xlsx --output results/test_part1.csv --seconds 60 --max-order 18 --seed 42
```

Résultat obtenu lors du run final :

```text
Conjectures réfutées : 87 / 100
Score total          : 1649.06 s
Résultats            : results/test_part1.csv
```

---

## Partie 2 — Architecture inspirée de FunSearch

Le fichier `graphbench_solver.py` correspond à la version améliorée.

Cette version reprend le moteur de recherche de la Partie 1 et ajoute :

- un score contextuel selon les invariants présents dans la conjecture ;
- des générateurs spécialisés ;
- des mutations enrichies ;
- un mécanisme de redémarrage lorsque la recherche stagne.

Commande d’exécution :

```bash
python graphbench_solver.py --input benchmark.xlsx --output results/test_part2.csv --seconds 60 --max-order 18 --seed 42
```

Résultat obtenu lors du run final :

```text
Conjectures réfutées : 94 / 100
Score total          : 867.76  s
Résultats            : results/test_part2.csv
```

---

## Résultats finaux

| Version | Conjectures réfutées | Taux de réussite | Score total |
|---|---:|---:|---:|
| Partie 1 — Heuristique simple | 87 / 100 | 87 % | 1649.06 s |
| Partie 2 — FunSearch-inspired | 94 / 100 | 94 % | 867.76  s |

La Partie 2 améliore donc la Partie 1 avec :

- 5 conjectures supplémentaires réfutées ;
-  une baisse du score total de 781.30 secondes ; ;
- une meilleure exploration grâce aux générateurs spécialisés, aux mutations enrichies et au mécanisme de redémarrage.

---

## Installation

Installer les dépendances Python avec :

```bash
pip install -r requirements.txt
```

Si besoin, les bibliothèques principales peuvent être installées directement avec :

```bash
pip install networkx numpy pandas openpyxl
```

---

## Exécution complète

Placez `benchmark.xlsx` dans le même dossier que les scripts Python.

### Lancer la Partie 1

```bash
python graphbench_part1.py --input benchmark.xlsx --output results/test_part1.csv --seconds 60 --max-order 18 
```

### Lancer la Partie 2

```bash
python graphbench_solver.py --input benchmark.xlsx --output results/test_part2.csv --seconds 60 --max-order 18 
```

---

## Format des fichiers CSV

Les fichiers CSV de sortie contiennent notamment les colonnes suivantes :

- `conjecture_id`
- `status`
- `validation_status`
- `score_variant`
- `counterexample_g6`
- `order`
- `size`
- `counter_x`
- `counter_y`
- `violation_margin`
- `cost`
- `execution_time_seconds`

Le champ `status` indique si une conjecture a été réfutée :

- `FOUND` : un contre-exemple a été trouvé ;
- `NOT_FOUND` : aucun contre-exemple n’a été trouvé dans la limite de temps.

Le champ `counterexample_g6` contient le graphe trouvé au format `graph6`.

---

## Scoring

Le score suit la règle suivante :

- si un contre-exemple est trouvé en `t` secondes, le coût vaut `t` ;
- si aucun contre-exemple n’est trouvé en 60 secondes, le coût vaut `120` ;
- le score total correspond à la somme des coûts sur toutes les conjectures.

L’objectif est donc de minimiser le score total.

---

## Exemple de sortie finale

Exemple de sortie obtenue pour la Partie 2 :

```text
============================================================
Conjectures refutees : 94 / 100
Score total          : 867.76 s
Resultats            : results/test_part2.csv
============================================================
```
---

## Auteurs

- Lydia AMROUCHE
- Hadi ISSA

Master 1 MIAGE — Optimisation Combinatoire  
Année universitaire 2025–2026
