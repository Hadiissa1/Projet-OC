#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GraphBench Challenge - Partie 2 : FunSearch-inspired amélioree

Objectif :
- Garder le moteur efficace de la Partie 1.
- Ajouter une couche Partie 2 utile :
  1. score principal stable,
  2. micro-bonus contextuel,
  3. candidats spécialisés selon les invariants,
  4. restart intelligent en cas de stagnation.

Auteur : Hadi ISSA - M1 MIAGE, Optimisation Combinatoire, 2025-2026
"""

import argparse
import ast
import csv
import math
import os
import random
import time
from fractions import Fraction
from itertools import combinations

import networkx as nx
import numpy as np
import pandas as pd


# ═══════════════════════════════════════════════════════════ ALIAS MAP ════════

ALIAS = {
    "order": "n",
    "size": "m",
    "diameter": "diam",
    "radius": "rad",
    "minimum_degree": "delta",
    "maximum_degree": "Delta",
    "average_degree": "avg",
    "triangle_number": "t",
    "clique_number": "omega",
    "domination_number": "gamma",
    "total_domination_number": "gamma_t",
    "independence_number": "alpha",
    "vertex_cover_number": "tau",
    "matching_number": "mu",
    "vertex_connectivity": "kappa",
    "edge_connectivity": "kappa_prime",
    "independent_domination_number": "iota",
    "density": "density",

    "n": "n",
    "m": "m",
    "diam": "diam",
    "rad": "rad",
    "delta": "delta",
    "Delta": "Delta",
    "avg": "avg",
    "t": "t",
    "omega": "omega",
    "gamma": "gamma",
    "gamma_t": "gamma_t",
    "alpha": "alpha",
    "tau": "tau",
    "mu": "mu",
    "iota": "iota",
    "kappa": "kappa",
    "kappa_prime": "kappa_prime",

    "randic_index": "randic_index",
    "harmonic_index": "harmonic_index",
    "first_zagreb_index": "first_zagreb_index",
    "second_zagreb_index": "second_zagreb_index",
    "proximity": "proximity",
    "remoteness": "remoteness",
    "largest_eigenvalue": "largest_eigenvalue",
    "largest_distance_eigenvalue": "largest_distance_eigenvalue",
    "second_smallest_laplace_eigenvalue": "second_smallest_laplace_eigenvalue",
}


def canon(name):
    return ALIAS.get(str(name).strip(), str(name).strip())


def get_inv(props, name):
    return props.get(canon(name), 0)


def involved_invariants(c):
    return {canon(c["x_name"]), canon(c["y_name"])}


# ═════════════════════════════════════════════════════════ INVARIANTS ═════════

def _alpha_backtrack(G):
    nodes = list(G.nodes())
    n = len(nodes)

    if n == 0:
        return 0

    idx = {v: i for i, v in enumerate(nodes)}
    adj = [frozenset(idx[u] for u in G.neighbors(v) if u in idx) for v in nodes]
    best = [0]

    def bt(candidates, size):
        if size + len(candidates) <= best[0]:
            return

        if not candidates:
            best[0] = max(best[0], size)
            return

        v = candidates[0]
        rest = candidates[1:]

        bt([u for u in rest if u not in adj[v]], size + 1)
        bt(rest, size)

    bt(list(range(n)), 0)
    return best[0]


def _alpha_greedy(G):
    remaining = set(G.nodes())
    count = 0

    while remaining:
        v = min(remaining, key=lambda x: G.degree(x))
        count += 1
        remaining -= {v} | set(G.neighbors(v))

    return count


def _gamma_exact(G):
    nodes = list(G.nodes())
    n = len(nodes)

    if n == 0:
        return 0

    idx = {v: i for i, v in enumerate(nodes)}
    closed = [
        frozenset({i} | {idx[u] for u in G.neighbors(v) if u in idx})
        for i, v in enumerate(nodes)
    ]
    all_n = frozenset(range(n))

    for size in range(1, n + 1):
        for S in combinations(range(n), size):
            if frozenset().union(*(closed[i] for i in S)) == all_n:
                return size

    return n


def _gamma_greedy(G):
    remaining = set(G.nodes())
    dominated = set()
    all_nodes = set(G.nodes())
    count = 0

    while dominated != all_nodes and remaining:
        best = max(remaining, key=lambda v: len(({v} | set(G.neighbors(v))) - dominated))
        dominated |= {best} | set(G.neighbors(best))
        remaining.discard(best)
        count += 1

    return count


def _gamma_t_exact(G):
    nodes = list(G.nodes())
    n = len(nodes)

    if n < 2:
        return n

    idx = {v: i for i, v in enumerate(nodes)}
    open_nbr = [
        frozenset(idx[u] for u in G.neighbors(v) if u in idx)
        for v in nodes
    ]
    all_n = frozenset(range(n))

    for size in range(2, n + 1):
        for S in combinations(range(n), size):
            if frozenset().union(*(open_nbr[i] for i in S)) == all_n:
                return size

    return n


def exact_parameters(G):
    if G.number_of_nodes() == 0:
        return {k: 0 for k in set(ALIAS.values())}

    n = G.number_of_nodes()
    m = G.number_of_edges()
    density = 2 * m / (n * (n - 1)) if n > 1 else 0.0

    props = {
        "n": n,
        "m": m,
        "density": density,
    }

    degrees = [d for _, d in G.degree()]
    props["delta"] = min(degrees) if degrees else 0
    props["Delta"] = max(degrees) if degrees else 0
    props["avg"] = 2 * m / n if n > 0 else 0

    connected = nx.is_connected(G)

    if connected and n > 1:
        try:
            ecc = nx.eccentricity(G)
            props["diam"] = max(ecc.values())
            props["rad"] = min(ecc.values())
        except Exception:
            props["diam"] = n
            props["rad"] = 0
    else:
        props["diam"] = n
        props["rad"] = 0

    try:
        tri = nx.triangles(G)
        props["t"] = sum(tri.values()) // 3
    except Exception:
        props["t"] = 0

    try:
        props["omega"] = max((len(c) for c in nx.find_cliques(G)), default=1)
    except Exception:
        props["omega"] = 1

    try:
        props["mu"] = len(nx.max_weight_matching(G, maxcardinality=True))
    except Exception:
        props["mu"] = m // max(n, 1)

    if n <= 21:
        props["alpha"] = _alpha_backtrack(G)
    else:
        props["alpha"] = _alpha_greedy(G)

    props["tau"] = n - props["alpha"]

    if n <= 16:
        props["gamma"] = _gamma_exact(G)
        props["gamma_t"] = _gamma_t_exact(G)
    else:
        props["gamma"] = _gamma_greedy(G)
        props["gamma_t"] = min(n, props["gamma"] + 1)

    # Approximation volontairement rapide.
    # Le calcul exact de iota ralentit fortement la recherche.
    props["iota"] = props["alpha"]

    if connected:
        try:
            props["kappa"] = nx.node_connectivity(G)
        except Exception:
            props["kappa"] = min(degrees) if degrees else 0

        try:
            props["kappa_prime"] = nx.edge_connectivity(G)
        except Exception:
            props["kappa_prime"] = min(degrees) if degrees else 0
    else:
        props["kappa"] = 0
        props["kappa_prime"] = 0

    if m > 0:
        props["randic_index"] = sum(
            1.0 / math.sqrt(max(G.degree(u), 1) * max(G.degree(v), 1))
            for u, v in G.edges()
        )
        props["harmonic_index"] = sum(
            2.0 / max(G.degree(u) + G.degree(v), 1)
            for u, v in G.edges()
        )
    else:
        props["randic_index"] = 0.0
        props["harmonic_index"] = 0.0

    props["first_zagreb_index"] = sum(d ** 2 for _, d in G.degree())
    props["second_zagreb_index"] = sum(G.degree(u) * G.degree(v) for u, v in G.edges())

    if connected and n > 1:
        try:
            dist = dict(nx.all_pairs_shortest_path_length(G))
            avg_d = [sum(dist[v].values()) / (n - 1) for v in G.nodes()]
            props["proximity"] = min(avg_d)
            props["remoteness"] = max(avg_d)
        except Exception:
            props["proximity"] = 0
            props["remoteness"] = n
    else:
        props["proximity"] = 0
        props["remoteness"] = n

    try:
        A = nx.adjacency_matrix(G).toarray().astype(float)
        ev_A = np.linalg.eigvalsh(A)
        props["largest_eigenvalue"] = float(np.max(np.abs(ev_A)))

        L = nx.laplacian_matrix(G).toarray().astype(float)
        ev_L = np.linalg.eigvalsh(L)
        props["second_smallest_laplace_eigenvalue"] = float(sorted(ev_L)[1]) if n > 1 else 0.0

        if connected:
            D = nx.floyd_warshall_numpy(G)
            ev_D = np.linalg.eigvalsh(D)
            props["largest_distance_eigenvalue"] = float(np.max(np.abs(ev_D)))
        else:
            props["largest_distance_eigenvalue"] = 0.0
    except Exception:
        props.setdefault("largest_eigenvalue", 0.0)
        props.setdefault("second_smallest_laplace_eigenvalue", 0.0)
        props.setdefault("largest_distance_eigenvalue", 0.0)

    return props


# ══════════════════════════════════════════════ CONJECTURE LOADING ════════════

def _parse_list(s):
    if isinstance(s, (list, tuple)):
        return list(s)

    try:
        return ast.literal_eval(str(s).strip())
    except Exception:
        return [str(s).strip()]


def _to_float(x):
    if isinstance(x, (int, float)):
        return float(x)

    s = str(x).strip()

    try:
        return float(Fraction(s))
    except Exception:
        try:
            return float(s)
        except Exception:
            return 1.0


def _parse_coeffs(s):
    if isinstance(s, (int, float)):
        return [_to_float(s)]

    if isinstance(s, (list, tuple, np.ndarray)):
        return [_to_float(x) for x in s]

    raw = str(s).strip()

    try:
        v = ast.literal_eval(raw)
        if isinstance(v, (int, float)):
            return [_to_float(v)]
        return [_to_float(x) for x in v]
    except Exception:
        try:
            return [float(Fraction(raw))]
        except Exception:
            return [1.0]


def load_conjectures(path):
    df = pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]

    # Colonnes attendues (au moins une variante pour chaque)
    required_columns = [
        ["Conjecture ID", "id", "ID"],
        ["Conjecture"],
        ["Subgroup", "subgroup"],
        ["X", "x"],
        ["Y", "y"],
        ["Sign", "sign"],
        ["Coefficients", "coefficients"],
        ["Intercept", "intercept"],
        ["Degree", "degree"]
    ]

    missing = []
    for variants in required_columns:
        if not any(col in df.columns for col in variants):
            missing.append(variants[0])
    if missing:
        raise ValueError(f"Le fichier benchmark est invalide. Colonnes manquantes : {', '.join(missing)}")

    conjectures = []

    for _, row in df.iterrows():
        try:
            cid = int(row.get("Conjecture ID", row.get("id", row.get("ID", 0))))
            text = str(row.get("Conjecture", ""))
            sg_raw = row.get("Subgroup", row.get("subgroup", "['connected']"))
            x_raw = str(row.get("X", row.get("x", "n"))).strip()
            y_raw = str(row.get("Y", row.get("y", "m"))).strip()
            sign = str(row.get("Sign", row.get("sign", "<="))).strip()
            coefs = _parse_coeffs(row.get("Coefficients", row.get("coefficients", [1.0])))
            intercept = _to_float(row.get("Intercept", row.get("intercept", 0.0)))
            degree = int(row.get("Degree", row.get("degree", 1)))
            subgroup = _parse_list(sg_raw)

            if not isinstance(subgroup, list):
                subgroup = [str(subgroup)]

            conjectures.append({
                "id": cid,
                "text": text,
                "subgroup": subgroup,
                "x_name": x_raw,
                "y_name": y_raw,
                "sign": sign,
                "coefficients": coefs,
                "intercept": intercept,
                "degree": degree,
            })

        except Exception as e:
            print(f"  Warning: row parse error: {e}")

    return conjectures


def eval_poly(c, x_val):
    coeffs = c["coefficients"]
    deg = c["degree"]
    result = float(c["intercept"])

    for i, coef in enumerate(coeffs):
        result += float(coef) * (float(x_val) ** (deg - i))

    return result


VIOLATION_THRESHOLD = 1e-9


def violation(c, props):
    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)

    if c["sign"] == "<=":
        return y_val - fx

    return fx - y_val


def _density(n, m):
    return 2 * m / (n * (n - 1)) if n > 1 else 0.0


# ═══════════════════════════════════════════ SCORE PRINCIPAL PARTIE 2 ════════

def score_base_v6(G, props, c):
    """
    Score principal stable, proche de la meilleure version empirique.
    """
    viol = violation(c, props)
    n = props.get("n", 1)

    if viol > VIOLATION_THRESHOLD:
        return 1e6 + 1e4 * viol - 0.01 * n

    m = props.get("m", 0)

    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)

    scale = 1 + abs(x_val) + abs(y_val) + abs(fx)
    margin_norm = viol / scale

    diam = props.get("diam", 0)
    alpha = props.get("alpha", 0)
    gamma = props.get("gamma", 0)
    omega = props.get("omega", 0)
    mu = props.get("mu", 0)
    kappa = props.get("kappa", 0)

    structure = (
        0.025 * diam
        + 0.020 * alpha
        + 0.015 * gamma
        + 0.010 * omega
        + 0.008 * mu
        + 0.005 * kappa
    )

    sg = (c["subgroup"][0] if c.get("subgroup") else "connected").lower()
    class_bonus = 0.0

    if sg == "tree":
        class_bonus = 0.020 * diam
    elif sg == "claw_free":
        class_bonus = 0.020 * omega
    elif sg == "bipartite":
        class_bonus = 0.015 * mu

    direction = 0.04 * viol
    signal = 0.03 * abs(viol) + 0.01 * abs(x_val)

    density = _density(n, m)
    complexity = 0.018 * n + 0.003 * m + 0.060 * abs(density - 0.45)

    return 700 * margin_norm + structure + class_bonus + direction + signal - complexity


def _contextual_bonus(G, props, c):
    """
    Micro-bonus contextuel.
    Il est volontairement faible pour ne pas casser l'ordre de sélection.
    """
    involved = involved_invariants(c)
    bonus = 0.0

    if "diam" in involved:
        bonus += 0.010 * props.get("diam", 0)

    if "rad" in involved:
        bonus += 0.010 * props.get("rad", 0)

    if "gamma" in involved:
        bonus += 0.010 * props.get("gamma", 0)

    if "gamma_t" in involved:
        bonus += 0.010 * props.get("gamma_t", 0)

    if "alpha" in involved:
        bonus += 0.010 * props.get("alpha", 0)

    if "tau" in involved:
        bonus += 0.008 * props.get("tau", 0)

    if "mu" in involved:
        bonus += 0.010 * props.get("mu", 0)

    if "omega" in involved:
        bonus += 0.010 * props.get("omega", 0)

    if "t" in involved:
        bonus += 0.0005 * props.get("t", 0)

    if "kappa" in involved:
        bonus += 0.010 * props.get("kappa", 0)

    if "kappa_prime" in involved:
        bonus += 0.010 * props.get("kappa_prime", 0)

    if "avg" in involved:
        bonus += 0.005 * props.get("avg", 0)

    if "Delta" in involved:
        bonus += 0.005 * props.get("Delta", 0)

    if "density" in involved:
        n = props.get("n", 1)
        m = props.get("m", 0)
        density = _density(n, m)
        bonus += 0.005 * (1.0 - abs(density - 0.45))

    if "remoteness" in involved:
        bonus += 0.005 * props.get("remoteness", 0)

    if "proximity" in involved:
        bonus -= 0.003 * props.get("proximity", 0)

    if "largest_eigenvalue" in involved:
        bonus += 0.005 * props.get("largest_eigenvalue", 0)

    if "largest_distance_eigenvalue" in involved:
        bonus += 0.003 * props.get("largest_distance_eigenvalue", 0)

    if "second_smallest_laplace_eigenvalue" in involved:
        bonus += 0.005 * props.get("second_smallest_laplace_eigenvalue", 0)

    return bonus


def choose_score(G, props, c):
    viol = violation(c, props)

    if viol > VIOLATION_THRESHOLD:
        return "score_base_v6_contextual", 1e6 + 1e4 * viol

    try:
        base_score = score_base_v6(G, props, c)
    except Exception:
        base_score = viol

    return "score_base_v6_contextual", base_score + _contextual_bonus(G, props, c)


# ════════════════════════════════════════════════════ GRAPH CLASS HELPERS ═════

def _is_claw_free(G):
    for v in G.nodes():
        nbrs = list(G.neighbors(v))

        if len(nbrs) >= 3:
            for a, b, c_n in combinations(nbrs, 3):
                if not (G.has_edge(a, b) or G.has_edge(b, c_n) or G.has_edge(a, c_n)):
                    return False

    return True


def is_in_class(G, subgroup):
    if G is None or G.number_of_nodes() == 0:
        return False

    sg = (subgroup[0] if subgroup else "connected").lower()

    if not nx.is_connected(G):
        return False

    if sg == "tree":
        return nx.is_tree(G)

    if sg == "claw_free":
        return _is_claw_free(G)

    if sg == "bipartite":
        return nx.is_bipartite(G)

    if sg == "planar":
        return nx.is_planar(G)

    return True


def repair_graph(G, subgroup, rng):
    if G.number_of_nodes() == 0:
        return nx.path_graph(3)

    G = G.copy()
    sg = (subgroup[0] if subgroup else "connected").lower()

    if not nx.is_connected(G):
        comps = list(nx.connected_components(G))
        for i in range(len(comps) - 1):
            u = rng.choice(list(comps[i]))
            v = rng.choice(list(comps[i + 1]))
            G.add_edge(u, v)

    if sg == "tree":
        if not nx.is_tree(G):
            G = nx.minimum_spanning_tree(G)

    elif sg == "claw_free":
        for _ in range(40):
            found_claw = False

            for v in list(G.nodes()):
                nbrs = list(G.neighbors(v))

                if len(nbrs) >= 3:
                    for a, b, c_n in combinations(nbrs, 3):
                        if not (G.has_edge(a, b) or G.has_edge(b, c_n) or G.has_edge(a, c_n)):
                            G.add_edge(a, b)
                            found_claw = True
                            break

                if found_claw:
                    break

            if not found_claw:
                break

    elif sg == "bipartite":
        if not nx.is_bipartite(G):
            color = {}
            root = list(G.nodes())[0]
            color[root] = 0
            queue = [root]

            while queue:
                v = queue.pop(0)

                for u in list(G.neighbors(v)):
                    if u not in color:
                        color[u] = 1 - color[v]
                        queue.append(u)
                    elif color[u] == color[v]:
                        G.remove_edge(v, u)

            if not nx.is_connected(G):
                comps = list(nx.connected_components(G))
                for i in range(len(comps) - 1):
                    u = list(comps[i])[0]
                    connected = False

                    for w in comps[i + 1]:
                        if color.get(w, 0) != color.get(u, 0):
                            G.add_edge(u, w)
                            connected = True
                            break

                    if not connected:
                        G.add_edge(u, list(comps[i + 1])[0])

    elif sg == "planar":
        if not nx.is_planar(G):
            G = nx.minimum_spanning_tree(G)

    return nx.convert_node_labels_to_integers(G)


# ═══════════════════════════════════════════ GENERATEURS SPECIALISES ══════════

def choose_mutation_profile(c):
    involved = involved_invariants(c)

    if "density" in involved:
        return "density"

    if "remoteness" in involved or "proximity" in involved:
        return "distance"

    if "diam" in involved or "rad" in involved:
        return "distance"

    if "largest_eigenvalue" in involved:
        return "spectral"

    if "largest_distance_eigenvalue" in involved:
        return "spectral"

    if "second_smallest_laplace_eigenvalue" in involved:
        return "spectral"

    if "randic_index" in involved or "harmonic_index" in involved:
        return "degree_index"

    if "first_zagreb_index" in involved or "second_zagreb_index" in involved:
        return "degree_index"

    if "gamma" in involved or "gamma_t" in involved:
        return "domination"

    return "general"


def lollipop_graph_custom(clique_size, path_len):
    G = nx.complete_graph(clique_size)
    last = clique_size - 1
    previous = last

    for i in range(path_len):
        new_v = clique_size + i
        G.add_edge(previous, new_v)
        previous = new_v

    return nx.convert_node_labels_to_integers(G)


def two_cliques_bridge(a, b, path_len):
    G1 = nx.complete_graph(a)
    G2 = nx.complete_graph(b)
    G = nx.disjoint_union(G1, G2)

    left = 0
    right = a

    previous = left
    next_id = a + b

    for _ in range(path_len):
        G.add_edge(previous, next_id)
        previous = next_id
        next_id += 1

    G.add_edge(previous, right)

    return nx.convert_node_labels_to_integers(G)


def generate_density_candidates(max_n, rng):
    graphs = []

    for n in range(5, max_n + 4):
        graphs.append(nx.path_graph(n))
        graphs.append(nx.cycle_graph(n))
        graphs.append(nx.star_graph(n - 1))
        graphs.append(nx.complete_graph(n))

        G = nx.complete_graph(n)
        edges = list(G.edges())
        rng.shuffle(edges)

        for e in edges[:max(1, n // 2)]:
            G.remove_edge(*e)
            if not nx.is_connected(G):
                G.add_edge(*e)

        graphs.append(G)

        for p in [0.12, 0.18, 0.25, 0.40, 0.60, 0.80]:
            H = nx.gnp_random_graph(n, p, seed=rng.randint(0, 999999))
            if nx.is_connected(H):
                graphs.append(H)

    return [
        nx.convert_node_labels_to_integers(G)
        for G in graphs
        if G.number_of_nodes() > 0 and nx.is_connected(G)
    ]


def generate_distance_candidates(max_n):
    graphs = []

    for n in range(5, max_n + 4):
        graphs.append(nx.path_graph(n))
        graphs.append(nx.star_graph(n - 1))
        graphs.append(nx.cycle_graph(n))

    for c in range(3, 9):
        for p in range(2, 12):
            G = lollipop_graph_custom(c, p)
            if G.number_of_nodes() <= max_n + 3:
                graphs.append(G)

    for a in range(3, 7):
        for b in range(3, 7):
            for p in range(1, 8):
                G = two_cliques_bridge(a, b, p)
                if G.number_of_nodes() <= max_n + 3:
                    graphs.append(G)

    return [
        nx.convert_node_labels_to_integers(G)
        for G in graphs
        if G.number_of_nodes() > 0 and nx.is_connected(G)
    ]


def generate_clawfree_candidates(max_n, rng):
    graphs = []

    for n in range(4, max_n + 1):
        for p in [0.20, 0.35, 0.50, 0.70]:
            base = nx.gnp_random_graph(n, p, seed=rng.randint(0, 999999))
            if base.number_of_edges() > 0:
                LG = nx.convert_node_labels_to_integers(nx.line_graph(base))
                if LG.number_of_nodes() > 0 and LG.number_of_nodes() <= max_n + 3 and nx.is_connected(LG):
                    graphs.append(LG)

    for n in range(4, max_n + 1):
        for base in [nx.path_graph(n), nx.cycle_graph(n), nx.complete_graph(min(n, 7))]:
            LG = nx.convert_node_labels_to_integers(nx.line_graph(base))
            if LG.number_of_nodes() > 0 and LG.number_of_nodes() <= max_n + 3 and nx.is_connected(LG):
                graphs.append(LG)

    for n in range(3, max_n + 1):
        graphs.append(nx.complete_graph(n))
        graphs.append(nx.cycle_graph(max(3, n)))

    return [
        nx.convert_node_labels_to_integers(G)
        for G in graphs
        if G.number_of_nodes() > 0 and nx.is_connected(G)
    ]


def generate_spectral_candidates(max_n, rng):
    graphs = []

    for n in range(5, max_n + 4):
        graphs.append(nx.complete_graph(n))
        graphs.append(nx.path_graph(n))
        graphs.append(nx.cycle_graph(n))
        graphs.append(nx.star_graph(n - 1))

        if n >= 4:
            graphs.append(nx.wheel_graph(n))

        for d in range(2, min(6, n)):
            if n * d % 2 == 0:
                try:
                    G = nx.random_regular_graph(d, n, seed=rng.randint(0, 999999))
                    if nx.is_connected(G):
                        graphs.append(G)
                except Exception:
                    pass

        for p in [0.20, 0.35, 0.50, 0.70, 0.85]:
            G = nx.gnp_random_graph(n, p, seed=rng.randint(0, 999999))
            if nx.is_connected(G):
                graphs.append(G)

    return [
        nx.convert_node_labels_to_integers(G)
        for G in graphs
        if G.number_of_nodes() > 0 and nx.is_connected(G)
    ]


def generate_degree_index_candidates(max_n, rng):
    graphs = []

    for n in range(5, max_n + 4):
        graphs.append(nx.star_graph(n - 1))
        graphs.append(nx.complete_graph(n))
        graphs.append(nx.path_graph(n))

        if n >= 4:
            graphs.append(nx.wheel_graph(n))

        for k in range(2, min(6, n)):
            G = nx.star_graph(n - 1)
            center = 0
            leaves = [v for v in G.nodes() if v != center]
            rng.shuffle(leaves)

            for i in range(min(k, len(leaves) - 1)):
                G.add_edge(leaves[i], leaves[i + 1])

            graphs.append(G)

        for p in [0.15, 0.30, 0.60, 0.85]:
            H = nx.gnp_random_graph(n, p, seed=rng.randint(0, 999999))
            if nx.is_connected(H):
                graphs.append(H)

    return [
        nx.convert_node_labels_to_integers(G)
        for G in graphs
        if G.number_of_nodes() > 0 and nx.is_connected(G)
    ]


def generate_domination_candidates(max_n, rng):
    graphs = []

    for n in range(5, max_n + 4):
        graphs.append(nx.path_graph(n))
        graphs.append(nx.cycle_graph(n))
        graphs.append(nx.star_graph(n - 1))

        # Couronne : cycle avec feuilles
        C = nx.cycle_graph(max(3, n // 2))
        base_nodes = list(C.nodes())
        next_id = C.number_of_nodes()

        for v in base_nodes:
            if next_id < max_n + 3:
                C.add_edge(v, next_id)
                next_id += 1

        graphs.append(nx.convert_node_labels_to_integers(C))

        # Plusieurs petites cliques reliées par un chemin
        if n >= 8:
            G = nx.Graph()
            current = 0
            anchors = []

            for _ in range(3):
                size = 3
                clique_nodes = list(range(current, current + size))
                for u, v in combinations(clique_nodes, 2):
                    G.add_edge(u, v)
                anchors.append(clique_nodes[0])
                current += size

            G.add_edge(anchors[0], anchors[1])
            G.add_edge(anchors[1], anchors[2])

            if G.number_of_nodes() <= max_n + 3:
                graphs.append(nx.convert_node_labels_to_integers(G))

    return [
        nx.convert_node_labels_to_integers(G)
        for G in graphs
        if G.number_of_nodes() > 0 and nx.is_connected(G)
    ]


def generate_specialized_candidates(c, sg, max_n, rng):
    profile = choose_mutation_profile(c)

    candidates = []

    if profile == "density":
        candidates.extend(generate_density_candidates(max_n, rng))

    elif profile == "distance":
        candidates.extend(generate_distance_candidates(max_n))

    elif profile == "spectral":
        candidates.extend(generate_spectral_candidates(max_n, rng))

    elif profile == "degree_index":
        candidates.extend(generate_degree_index_candidates(max_n, rng))

    elif profile == "domination":
        candidates.extend(generate_domination_candidates(max_n, rng))

    if sg == "claw_free":
        candidates.extend(generate_clawfree_candidates(max_n, rng))

    rng.shuffle(candidates)
    return candidates


# ════════════════════════════════════════════════════ INITIAL BANK ════════════

def _class_specific(sg, max_n, rng):
    graphs = []

    if sg == "tree":
        for n in range(2, min(max_n + 1, 18)):
            graphs.append(nx.path_graph(n))
            graphs.append(nx.star_graph(n))

            T = nx.path_graph(n)
            for v in range(n):
                if rng.random() < 0.4 and T.number_of_nodes() < max_n:
                    T.add_edge(v, T.number_of_nodes())
            graphs.append(T.copy())

        for n in range(4, min(max_n + 1, 15)):
            for _ in range(8):
                seq = [rng.randint(0, n - 1) for _ in range(n - 2)]
                try:
                    graphs.append(nx.from_prufer_sequence(seq))
                except Exception:
                    pass

    elif sg == "bipartite":
        for a in range(1, 9):
            for b in range(a, 9):
                if a + b <= max_n + 2:
                    graphs.append(nx.complete_bipartite_graph(a, b))

        for n in range(4, min(max_n + 1, 15)):
            half = n // 2
            G = nx.bipartite.random_graph(half, n - half, 0.5, seed=rng.randint(0, 99999))
            if nx.is_connected(G):
                graphs.append(nx.convert_node_labels_to_integers(G))

    elif sg == "claw_free":
        graphs.extend(generate_clawfree_candidates(max_n, rng))

    elif sg == "planar":
        for n in range(2, min(max_n + 1, 16)):
            graphs.append(nx.path_graph(n))

            if n >= 3:
                graphs.append(nx.cycle_graph(n))

            graphs.append(nx.ladder_graph(n))

        for k in range(2, 6):
            graphs.append(nx.grid_2d_graph(k, k))

        graphs += [nx.dodecahedral_graph(), nx.octahedral_graph(), nx.cubical_graph()]

    else:
        for n in range(2, min(max_n + 1, 12)):
            graphs.append(nx.path_graph(n))

            if n >= 3:
                graphs.append(nx.cycle_graph(n))

            graphs.append(nx.complete_graph(n))
            graphs.append(nx.star_graph(n))

            if n >= 3:
                graphs.append(nx.wheel_graph(n))

    result = []

    for G in graphs:
        G = nx.convert_node_labels_to_integers(G)

        if G.number_of_nodes() >= 1 and nx.is_connected(G):
            result.append(G.copy())

    return result


def generate_initial_bank(max_n, rng):
    bank = []

    try:
        for G in nx.graph_atlas_g():
            if 1 <= G.number_of_nodes() <= 7 and nx.is_connected(G):
                bank.append(G.copy())
    except Exception:
        pass

    for n in range(2, min(max_n + 1, 19)):
        bank.append(nx.path_graph(n))

        if n >= 3:
            bank.append(nx.cycle_graph(n))

        bank.append(nx.complete_graph(n))
        bank.append(nx.star_graph(n))

        if n >= 3:
            bank.append(nx.wheel_graph(n))

        bank.append(nx.ladder_graph(n))

    for a in range(1, 7):
        for b in range(a, 7):
            if a + b <= max_n + 2:
                bank.append(nx.complete_bipartite_graph(a, b))

    for n in range(4, min(max_n + 1, 14)):
        for _ in range(10):
            seq = [rng.randint(0, n - 1) for _ in range(n - 2)]

            try:
                bank.append(nx.from_prufer_sequence(seq))
            except Exception:
                pass

    for n in range(4, min(max_n + 1, 14)):
        for p in [0.10, 0.18, 0.25, 0.35, 0.50, 0.70, 0.85]:
            for _ in range(5):
                G = nx.gnp_random_graph(n, p, seed=rng.randint(0, 99999))
                if nx.is_connected(G):
                    bank.append(G.copy())

    for n in range(4, min(max_n + 1, 13)):
        for d in range(2, min(6, n)):
            if n * d % 2 == 0:
                try:
                    G = nx.random_regular_graph(d, n, seed=rng.randint(0, 99999))
                    if nx.is_connected(G):
                        bank.append(G.copy())
                except Exception:
                    pass

    for n in range(3, 9):
        for base in [nx.path_graph(n), nx.cycle_graph(n), nx.complete_graph(n)]:
            LG = nx.convert_node_labels_to_integers(nx.line_graph(base))
            if nx.is_connected(LG):
                bank.append(LG)

    for G in [
        nx.petersen_graph(),
        nx.dodecahedral_graph(),
        nx.icosahedral_graph(),
        nx.octahedral_graph(),
        nx.cubical_graph(),
    ]:
        bank.append(G)

    result = []

    for G in bank:
        G = nx.convert_node_labels_to_integers(G)

        if 1 <= G.number_of_nodes() <= max_n + 3 and nx.is_connected(G):
            result.append(G.copy())

    return result


# ════════════════════════════════════════════════════ MUTATIONS ENRICHIES ═════

def _mutate_tree(G, max_n, rng):
    G = G.copy()
    nodes = list(G.nodes())
    n = len(nodes)

    ops = ["add_leaf", "subdivide", "move_leaf", "add_leaf"]

    if n > 3:
        ops.append("remove_leaf")

    op = rng.choice(ops)

    if op == "add_leaf" and n < max_n + 3:
        G.add_edge(rng.choice(nodes), max(nodes) + 1)

    elif op == "remove_leaf":
        leaves = [v for v in nodes if G.degree(v) == 1]
        if leaves:
            G.remove_node(rng.choice(leaves))

    elif op == "subdivide":
        edges = list(G.edges())
        if edges and n < max_n + 3:
            u, v = rng.choice(edges)
            new_v = max(G.nodes()) + 1
            G.remove_edge(u, v)
            G.add_edge(u, new_v)
            G.add_edge(new_v, v)

    elif op == "move_leaf":
        leaves = [v for v in nodes if G.degree(v) == 1]
        non_leaves = [v for v in nodes if G.degree(v) > 1]

        if leaves and non_leaves:
            leaf = rng.choice(leaves)
            candidates = [v for v in non_leaves if v != leaf]

            if candidates:
                new_p = rng.choice(candidates)
                old_p = list(G.neighbors(leaf))[0]

                if new_p != old_p:
                    G.remove_edge(leaf, old_p)
                    G.add_edge(leaf, new_p)

    return G


def mutate(G, subgroup, max_n, rng, steps=1, profile="general"):
    sg = (subgroup[0] if subgroup else "connected").lower()

    for _ in range(steps):
        G = _single_mutate(G, sg, max_n, rng, profile)

    return G


def _single_mutate(G, sg, max_n, rng, profile="general"):
    G = G.copy()
    n = G.number_of_nodes()
    nodes = list(G.nodes())
    edges = list(G.edges())

    if sg == "tree":
        return _mutate_tree(G, max_n, rng)

    ops = [
        "add_edge",
        "remove_edge",
        "toggle_edge",
        "add_vertex",
        "subdivide",
        "add_path",
        "add_clique",
        "add_twins",
        "densify_local",
        "sparsify_local",
        "attach_star",
        "attach_cycle",
        "double_subdivide",
    ]

    if profile == "density":
        ops += ["add_edge", "add_edge", "densify_local", "sparsify_local"]

    elif profile == "distance":
        ops += ["add_path", "subdivide", "double_subdivide", "attach_star"]

    elif profile == "spectral":
        ops += ["add_edge", "densify_local", "add_clique", "line_graph"]

    elif profile == "degree_index":
        ops += ["attach_star", "add_twins", "add_clique", "densify_local"]

    elif profile == "domination":
        ops += ["add_path", "attach_star", "add_twins", "sparsify_local"]

    if n > 2:
        ops.append("remove_vertex")

    if n <= 10:
        ops.append("line_graph")

    op = rng.choice(ops)

    if op == "add_edge":
        ne = list(nx.non_edges(G))
        if ne:
            G.add_edge(*rng.choice(ne))

    elif op == "remove_edge" and edges:
        rng.shuffle(edges)

        for u, v in edges:
            G.remove_edge(u, v)

            if nx.is_connected(G):
                break

            G.add_edge(u, v)

    elif op == "toggle_edge":
        ne = list(nx.non_edges(G))

        if ne and (not edges or rng.random() < 0.5):
            G.add_edge(*rng.choice(ne))

        elif edges:
            u, v = rng.choice(edges)
            G.remove_edge(u, v)

            if not nx.is_connected(G):
                G.add_edge(u, v)

    elif op == "add_vertex" and n < max_n + 3:
        new_v = max(nodes) + 1
        k = rng.randint(1, min(3, n))

        for nb in rng.sample(nodes, k):
            G.add_edge(new_v, nb)

    elif op == "remove_vertex" and n > 2:
        rng.shuffle(nodes)

        for v in nodes:
            H = G.copy()
            H.remove_node(v)

            if H.number_of_nodes() > 0 and nx.is_connected(H):
                G = H
                break

    elif op == "subdivide" and edges and n < max_n + 3:
        u, v = rng.choice(edges)
        new_v = max(G.nodes()) + 1
        G.remove_edge(u, v)
        G.add_edge(u, new_v)
        G.add_edge(new_v, v)

    elif op == "double_subdivide" and edges and n < max_n + 2:
        u, v = rng.choice(edges)
        nv1 = max(G.nodes()) + 1
        nv2 = nv1 + 1
        G.remove_edge(u, v)
        G.add_edge(u, nv1)
        G.add_edge(nv1, nv2)
        G.add_edge(nv2, v)

    elif op == "add_path" and n < max_n + 3:
        attach = rng.choice(nodes)
        new_v = max(G.nodes()) + 1
        G.add_edge(attach, new_v)

        if n < max_n + 2 and rng.random() < 0.5:
            G.add_edge(new_v, new_v + 1)

    elif op == "add_clique" and n < max_n + 2:
        attach = rng.choice(nodes)
        nv1 = max(G.nodes()) + 1
        nv2 = nv1 + 1
        G.add_edge(attach, nv1)
        G.add_edge(attach, nv2)
        G.add_edge(nv1, nv2)

    elif op == "add_twins" and n < max_n + 3:
        v = rng.choice(nodes)
        nbrs = list(G.neighbors(v))
        new_v = max(G.nodes()) + 1

        for u in nbrs:
            G.add_edge(new_v, u)

        if rng.random() < 0.5:
            G.add_edge(new_v, v)

        if new_v not in G.nodes() or G.degree(new_v) == 0:
            G.add_edge(new_v, rng.choice(nodes))

    elif op == "densify_local":
        subset = rng.sample(nodes, min(5, n))

        for i in range(len(subset)):
            for j in range(i + 1, len(subset)):
                if rng.random() < 0.4 and not G.has_edge(subset[i], subset[j]):
                    G.add_edge(subset[i], subset[j])

    elif op == "sparsify_local":
        removable = list(edges)
        rng.shuffle(removable)
        removed = 0

        for u, v in removable:
            if removed >= 2:
                break

            G.remove_edge(u, v)

            if nx.is_connected(G):
                removed += 1
            else:
                G.add_edge(u, v)

    elif op == "attach_star" and n < max_n + 2:
        attach = rng.choice(nodes)
        nv1 = max(G.nodes()) + 1
        nv2 = nv1 + 1
        G.add_edge(attach, nv1)
        G.add_edge(attach, nv2)

    elif op == "attach_cycle" and n <= max_n - 1:
        attach = rng.choice(nodes)
        base = max(G.nodes()) + 1

        for i in range(4):
            G.add_edge(base + i, base + (i + 1) % 4)

        G.add_edge(attach, base)

    elif op == "line_graph" and n <= 10:
        LG = nx.convert_node_labels_to_integers(nx.line_graph(G))

        if nx.is_connected(LG) and 1 <= LG.number_of_nodes() <= max_n + 3:
            G = LG

    return G


def _fresh_candidate(sg, max_n, rng):
    n = rng.randint(3, max(4, max_n - 2))

    if sg == "tree":
        if rng.random() < 0.5:
            seq = [rng.randint(0, n - 1) for _ in range(max(1, n - 2))]

            try:
                return nx.from_prufer_sequence(seq)
            except Exception:
                return nx.path_graph(n)

        return nx.star_graph(n)

    if sg == "bipartite":
        a = rng.randint(1, max(2, n // 2))
        b = n - a
        G = nx.bipartite.random_graph(a, b, rng.random() * 0.7 + 0.1, seed=rng.randint(0, 99999))

        if nx.is_connected(G):
            return nx.convert_node_labels_to_integers(G)

        return nx.complete_bipartite_graph(a, b)

    if sg == "claw_free":
        base = nx.gnm_random_graph(n, rng.randint(n, 2 * n), seed=rng.randint(0, 99999))
        LG = nx.convert_node_labels_to_integers(nx.line_graph(base))

        if LG.number_of_nodes() > 0 and nx.is_connected(LG):
            return LG

        return nx.complete_graph(n)

    if sg == "planar":
        choices = [
            nx.path_graph(n),
            nx.cycle_graph(n),
            nx.ladder_graph(max(2, n // 2)),
        ]
        return nx.convert_node_labels_to_integers(rng.choice(choices))

    p = rng.random() * 0.5 + 0.2
    G = nx.gnp_random_graph(n, p, seed=rng.randint(0, 99999))

    if nx.is_connected(G):
        return nx.convert_node_labels_to_integers(G)

    return nx.path_graph(n)


# ══════════════════════════════════════════════════════ MAIN SOLVER ═══════════

def _make_result(c, G, props, viol, cost, status, score_variant, exec_time=None):
    g6 = ""
    order = 0
    size = 0
    x_val = 0.0
    y_val = 0.0

    if G is not None and G.number_of_nodes() > 0:
        try:
            g6 = nx.to_graph6_bytes(G, header=False).decode().strip()
        except Exception:
            g6 = ""

        order = G.number_of_nodes()
        size = G.number_of_edges()

    if props is not None:
        x_val = get_inv(props, c["x_name"])
        y_val = get_inv(props, c["y_name"])

    v_margin = round(viol, 6) if viol is not None and not math.isnan(viol) else -999

    return {
        "conjecture_id": c["id"],
        "status": status,
        "validation_status": "VALID_COUNTEREXAMPLE" if status == "FOUND" else "NOT_FOUND",
        "score_variant": score_variant,
        "counterexample_g6": g6,
        "order": order,
        "size": size,
        "counter_x": round(float(x_val), 6),
        "counter_y": round(float(y_val), 6),
        "violation_margin": v_margin,
        "cost": round(float(cost), 4),
        "execution_time_seconds": round(float(exec_time if exec_time is not None else cost), 4),
    }


def _valid_counterexample(G, props, c, subgroup):
    if G is None or props is None:
        return False

    if not is_in_class(G, subgroup):
        return False

    return violation(c, props) > VIOLATION_THRESHOLD


def add_candidate_to_population(G, c, subgroup, rng, population, max_n, t0):
    try:
        G = repair_graph(G, subgroup, rng)

        if G.number_of_nodes() < 1 or G.number_of_nodes() > max_n + 3:
            return None

        if not is_in_class(G, subgroup):
            return None

        props = exact_parameters(G)
        viol = violation(c, props)

        if viol > VIOLATION_THRESHOLD:
            return ("FOUND", G, props, viol)

        name, s = choose_score(G, props, c)
        population.append((s, G.copy(), props, name))

        return ("ADDED", G, props, viol)

    except Exception:
        return None


def solve_one(c, initial_bank, seconds, max_n, rng):
    t0 = time.time()
    subgroup = c["subgroup"]
    sg = (subgroup[0] if subgroup else "connected").lower()
    profile = choose_mutation_profile(c)

    class_bank = [G for G in initial_bank if is_in_class(G, subgroup)]

    if len(class_bank) < 5:
        class_bank = _class_specific(sg, max_n, rng)

    specialized_start = generate_specialized_candidates(c, sg, max_n, rng)
    class_bank = (class_bank + specialized_start[:80])[:420]

    best_score = -math.inf
    best_name = "score_base_v6_contextual"
    best_G = None
    best_props = None
    best_viol = -math.inf

    population = []

    # Phase 1 : banque initiale + candidats spécialisés légers.
    for G in class_bank:
        if time.time() - t0 > seconds * 0.30:
            break

        try:
            G = repair_graph(G, subgroup, rng)

            if G.number_of_nodes() < 1 or G.number_of_nodes() > max_n + 3:
                continue

            if not is_in_class(G, subgroup):
                continue

            props = exact_parameters(G)
            viol = violation(c, props)

            if _valid_counterexample(G, props, c, subgroup):
                return _make_result(c, G, props, viol, time.time() - t0, "FOUND", "initial_or_specialized_bank")

            name, s = choose_score(G, props, c)
            population.append((s, G.copy(), props, name))

            if s > best_score:
                best_score = s
                best_G = G.copy()
                best_props = props
                best_viol = viol
                best_name = name

        except Exception:
            pass

    population.sort(key=lambda x: x[0], reverse=True)
    population = population[:140]

    TOP_BRACKETS = [12, 20, 35, 60]
    last_improvement = time.time()
    last_restart = time.time()

    # Phase 2 : recherche évolutive + restart intelligent.
    while time.time() - t0 < seconds:
        now = time.time()
        elapsed = now - t0

        # Restart intelligent si stagnation.
        if now - last_improvement > 8.0 and now - last_restart > 6.0:
            population = population[:20]

            extra = generate_specialized_candidates(c, sg, max_n, rng)

            if not extra:
                extra = [_fresh_candidate(sg, max_n, rng) for _ in range(60)]

            for G0 in extra[:90]:
                if time.time() - t0 >= seconds:
                    break

                result = add_candidate_to_population(G0, c, subgroup, rng, population, max_n, t0)

                if result and result[0] == "FOUND":
                    _, GG, pp, vv = result
                    return _make_result(
                        c,
                        GG,
                        pp,
                        vv,
                        time.time() - t0,
                        "FOUND",
                        "restart_specialized",
                    )

            population.sort(key=lambda x: x[0], reverse=True)
            population = population[:140]

            last_restart = time.time()
            last_improvement = time.time()

        r = rng.random()

        # En fin de temps, augmenter l'exploitation des meilleurs.
        if elapsed > seconds * 0.70:
            exploit_prob = 0.88
        else:
            exploit_prob = 0.78

        if r < exploit_prob and population:
            k = min(len(population), rng.choice(TOP_BRACKETS))
            idx = int(rng.random() ** 1.5 * k)
            parent = population[idx][1]

        elif r < 0.95:
            try:
                # En phase tardive, créer un candidat spécialisé.
                if elapsed > seconds * 0.55 and rng.random() < 0.70:
                    extra = generate_specialized_candidates(c, sg, max_n, rng)
                    parent = rng.choice(extra) if extra else _fresh_candidate(sg, max_n, rng)
                else:
                    parent = _fresh_candidate(sg, max_n, rng)

                parent = nx.convert_node_labels_to_integers(parent)

            except Exception:
                parent = rng.choice(class_bank) if class_bank else nx.path_graph(4)

        else:
            parent = rng.choice(class_bank) if class_bank else nx.path_graph(4)

        r2 = rng.random()

        if elapsed > seconds * 0.70:
            steps = 1 if r2 < 0.55 else (2 if r2 < 0.88 else 3)
        else:
            steps = 1 if r2 < 0.70 else (2 if r2 < 0.90 else 3)

        try:
            H = mutate(parent, subgroup, max_n, rng, steps, profile=profile)
            H = repair_graph(H, subgroup, rng)

            if H.number_of_nodes() < 1 or H.number_of_nodes() > max_n + 3:
                continue

            if not is_in_class(H, subgroup):
                continue

            props = exact_parameters(H)
            viol = violation(c, props)

            if _valid_counterexample(H, props, c, subgroup):
                name, _ = choose_score(H, props, c)
                return _make_result(c, H, props, viol, time.time() - t0, "FOUND", name)

            fn_name, s = choose_score(H, props, c)

            if len(population) < 140 or s > population[-1][0]:
                population.append((s, H.copy(), props, fn_name))
                population.sort(key=lambda x: x[0], reverse=True)
                population = population[:140]

            if s > best_score:
                best_score = s
                best_G = H.copy()
                best_props = props
                best_viol = viol
                best_name = fn_name
                last_improvement = time.time()

        except Exception:
            pass

    exec_time = time.time() - t0

    return _make_result(
        c,
        best_G,
        best_props,
        best_viol if best_viol != -math.inf else -999,
        120,
        "NOT_FOUND",
        best_name,
        exec_time,
    )


def solve(input_path, output_path, seconds, max_order, seed=42):
    rng = random.Random(seed)

    print(f"Loading conjectures from {input_path} ...")
    conjectures = load_conjectures(input_path)
    print(f"  {len(conjectures)} conjectures loaded.")

    print(f"Generating initial bank (max_order={max_order}) ...")
    bank = generate_initial_bank(max_order, rng)
    print(f"  {len(bank)} graphs in bank.")

    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fieldnames = [
        "conjecture_id",
        "status",
        "validation_status",
        "score_variant",
        "counterexample_g6",
        "order",
        "size",
        "counter_x",
        "counter_y",
        "violation_margin",
        "cost",
        "execution_time_seconds",
    ]

    total_cost = 0.0
    found = 0

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i, c in enumerate(conjectures):
            label = f"[{i + 1}/{len(conjectures)}] #{c['id']}"
            print(f"{label:20s} {c['text'][:55]}", end=" ... ", flush=True)

            try:
                res = solve_one(c, bank, seconds, max_order, random.Random(seed + c["id"]))
            except Exception as e:
                print(f"ERROR: {e}")
                res = _make_result(c, None, None, -999, 120, "NOT_FOUND", "error")

            total_cost += res["cost"]

            if res["status"] == "FOUND":
                found += 1
                print(
                    f"FOUND  t={res['cost']:.2f}s  "
                    f"variant={res['score_variant']}  "
                    f"margin={res['violation_margin']}"
                )
            else:
                print("NOT FOUND  (cost=120)")

            writer.writerow(res)
            f.flush()

    print(f"\n{'=' * 60}")
    print(f"Conjectures refutees : {found} / {len(conjectures)}")
    print(f"Score total          : {total_cost:.2f} s")
    print(f"Resultats            : {output_path}")
    print(f"{'=' * 60}")


def main():
    p = argparse.ArgumentParser(
        description="GraphBench Partie 2 - FunSearch-inspired avec strategies specialisees"
    )

    p.add_argument("--input", required=True, help="benchmark.xlsx")
    p.add_argument("--output", required=True, help="fichier CSV resultats")
    p.add_argument("--seconds", type=float, default=60.0, help="limite temps / conjecture")
    p.add_argument("--max-order", type=int, default=18, help="ordre max des graphes")
    p.add_argument("--seed", type=int, default=42, help="graine aleatoire")

    args = p.parse_args()
    solve(args.input, args.output, args.seconds, args.max_order, args.seed)


if __name__ == "__main__":
    main()