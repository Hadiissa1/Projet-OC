#!/usr/bin/env python3
"""
GraphBench Challenge - Partie 2 : Architecture FunSearch-inspired
Portfolio de 9 fonctions de score evaluees en parallele + selection automatique
Mutations enrichies + population evolutive de 140 candidats
Auteur : Hadi ISSA  M1 MIAGE, Optimisation Combinatoire, 2025-2026
"""

import argparse, ast, csv, math, os, random, time
from fractions import Fraction
from itertools import combinations
import networkx as nx
import numpy as np
import pandas as pd


# ═══════════════════════════════════════════════════════════ ALIAS MAP ════════

ALIAS = {
    "order": "n", "size": "m", "diameter": "diam", "radius": "rad",
    "minimum_degree": "delta", "maximum_degree": "Delta",
    "average_degree": "avg", "triangle_number": "t",
    "clique_number": "omega", "domination_number": "gamma",
    "total_domination_number": "gamma_t", "independence_number": "alpha",
    "vertex_cover_number": "tau", "matching_number": "mu",
    "vertex_connectivity": "kappa", "edge_connectivity": "kappa_prime",
    "independent_domination_number": "alpha",
    "density": "density",
    "n": "n", "m": "m", "diam": "diam", "rad": "rad",
    "delta": "delta", "Delta": "Delta", "avg": "avg", "t": "t",
    "omega": "omega", "gamma": "gamma", "gamma_t": "gamma_t",
    "alpha": "alpha", "tau": "tau", "mu": "mu",
    "kappa": "kappa", "kappa_prime": "kappa_prime",
    "randic_index": "randic_index", "harmonic_index": "harmonic_index",
    "first_zagreb_index": "first_zagreb_index",
    "second_zagreb_index": "second_zagreb_index",
    "proximity": "proximity", "remoteness": "remoteness",
    "largest_eigenvalue": "largest_eigenvalue",
    "largest_distance_eigenvalue": "largest_distance_eigenvalue",
    "second_smallest_laplace_eigenvalue": "second_smallest_laplace_eigenvalue",
}


def get_inv(props, name):
    return props.get(ALIAS.get(name.strip(), name.strip()), 0)


# ═════════════════════════════════════════════════════════ INVARIANTS ═════════

def _alpha_backtrack(G):
    nodes = list(G.nodes())
    n = len(nodes)
    if n == 0:
        return 0
    idx = {v: i for i, v in enumerate(nodes)}
    adj = [frozenset(idx[u] for u in G.neighbors(v) if u in idx) for v in nodes]
    best = [1]

    def bt(candidates, size):
        if size + len(candidates) <= best[0]:
            return
        if not candidates:
            if size > best[0]:
                best[0] = size
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
    closed = [frozenset({i} | {idx[u] for u in G.neighbors(v) if u in idx})
              for i, v in enumerate(nodes)]
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
    open_nbr = [frozenset(idx[u] for u in G.neighbors(v) if u in idx) for v in nodes]
    all_n = frozenset(range(n))
    for size in range(2, n + 1):
        for S in combinations(range(n), size):
            if frozenset().union(*(open_nbr[i] for i in S)) == all_n:
                return size
    return n


def exact_parameters(G):
    if G.number_of_nodes() == 0:
        return {k: 0 for k in ALIAS.values()}

    n = G.number_of_nodes()
    m = G.number_of_edges()
    density = 2 * m / (n * (n - 1)) if n > 1 else 0.0
    props = {"n": n, "m": m, "density": density}

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
        matching = nx.max_weight_matching(G, maxcardinality=True)
        props["mu"] = len(matching)
    except Exception:
        props["mu"] = m // max(n, 1)

    if n <= 20:
        props["alpha"] = _alpha_backtrack(G)
    else:
        props["alpha"] = _alpha_greedy(G)
    props["tau"] = n - props["alpha"]

    if n <= 15:
        props["gamma"] = _gamma_exact(G)
        props["gamma_t"] = _gamma_t_exact(G)
    else:
        props["gamma"] = _gamma_greedy(G)
        props["gamma_t"] = props["gamma"] + 1

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
            for u, v in G.edges())
        props["harmonic_index"] = sum(
            2.0 / max(G.degree(u) + G.degree(v), 1)
            for u, v in G.edges())
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
    """Convert int, float, or rational string like '3/4' to float."""
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
                "id": cid, "text": text, "subgroup": subgroup,
                "x_name": x_raw, "y_name": y_raw, "sign": sign,
                "coefficients": coefs, "intercept": intercept, "degree": degree,
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
    return (y_val - fx) if c["sign"] == "<=" else (fx - y_val)


# ═══════════════════════════════════════════ PORTFOLIO DE 9 FONCTIONS ════════
# Inspiration FunSearch : chaque fonction est une strategie independante.
# La selection automatique prend le max des 9 a chaque evaluation.

def _density(n, m):
    return 2 * m / (n * (n - 1)) if n > 1 else 0


def _v0_violation_pure(G, props, c):
    """v0: violation brute uniquement."""
    return violation(c, props)


def _v1_normalized_margin(G, props, c):
    """v1: marge normalisee x900 pour amplifier le signal."""
    viol = violation(c, props)
    if viol > VIOLATION_THRESHOLD:
        return 1e6 + 900 * viol
    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)
    scale = 1 + abs(x_val) + abs(y_val) + abs(fx)
    return 900 * viol / scale


def _v2_structure_guided(G, props, c):
    """v2: marge + bonus structure (diametre, gamma, alpha)."""
    viol = violation(c, props)
    n = props.get("n", 1)
    if viol > VIOLATION_THRESHOLD:
        return 1e6 + 500 * viol - 0.01 * n
    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)
    scale = 1 + abs(x_val) + abs(y_val) + abs(fx)
    m_norm = viol / scale
    diam = props.get("diam", 0)
    alpha = props.get("alpha", 0)
    gamma = props.get("gamma", 0)
    omega = props.get("omega", 0)
    mu = props.get("mu", 0)
    bonus = 0.03 * diam + 0.02 * alpha + 0.015 * gamma + 0.01 * omega + 0.005 * mu
    return 600 * m_norm + bonus - 0.015 * n


def _v3_conjecture_direction(G, props, c):
    """v3: marge + direction selon signe de la conjecture."""
    viol = violation(c, props)
    n = props.get("n", 1)
    if viol > VIOLATION_THRESHOLD:
        return 1e6 + 700 * viol - 0.01 * n
    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)
    scale = 1 + abs(x_val) + abs(y_val) + abs(fx)
    m_norm = viol / scale
    # penalise si on est tres loin de la violation
    direction_bonus = 0.05 * viol  # viol is negative here, guides toward 0
    return 550 * m_norm + direction_bonus


def _v4_class_specialized(G, props, c):
    """v4: bonus specifique a la classe de graphes."""
    viol = violation(c, props)
    n = props.get("n", 1)
    m = props.get("m", 0)
    if viol > VIOLATION_THRESHOLD:
        return 1e6 + 600 * viol - 0.01 * n
    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)
    scale = 1 + abs(x_val) + abs(y_val) + abs(fx)
    m_norm = viol / scale

    sg = (c["subgroup"][0] if c.get("subgroup") else "connected").lower()
    bonus = 0.0
    if sg == "tree":
        # arbres longs avec beaucoup de feuilles
        diam = props.get("diam", 0)
        leaves = sum(1 for v in G.nodes() if G.degree(v) == 1)
        bonus = 0.04 * diam + 0.02 * leaves
    elif sg == "claw_free":
        # graphes denses, cliques larges
        omega = props.get("omega", 0)
        bonus = 0.05 * omega + 0.02 * m / max(n, 1)
    elif sg == "bipartite":
        # bipartis complets ou denses
        mu = props.get("mu", 0)
        bonus = 0.04 * mu
    elif sg == "planar":
        # planaires avec grand diametre
        diam = props.get("diam", 0)
        bonus = 0.03 * diam
    else:
        # connexes generaux
        alpha = props.get("alpha", 0)
        gamma = props.get("gamma", 0)
        bonus = 0.02 * alpha + 0.015 * gamma

    density = _density(n, m)
    return 550 * m_norm + bonus - 0.015 * n - 0.05 * abs(density - 0.4)


def _v5_hard_invariant_guided(G, props, c):
    """v5: bonus si invariants difficiles impliques (alpha, gamma, kappa)."""
    viol = violation(c, props)
    n = props.get("n", 1)
    if viol > VIOLATION_THRESHOLD:
        return 1e6 + 800 * viol - 0.005 * n

    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)
    scale = 1 + abs(x_val) + abs(y_val) + abs(fx)
    m_norm = viol / scale

    # Bonus for invariants that are hard and relevant
    alpha = props.get("alpha", 0)
    gamma = props.get("gamma", 0)
    tau = props.get("tau", 0)
    kappa = props.get("kappa", 0)
    mu = props.get("mu", 0)
    omega = props.get("omega", 0)

    hard_names = {"alpha", "gamma", "gamma_t", "tau", "kappa", "kappa_prime",
                  "independence_number", "domination_number", "vertex_connectivity"}
    in_hard = (ALIAS.get(c["x_name"], c["x_name"]) in hard_names or
               ALIAS.get(c["y_name"], c["y_name"]) in hard_names)

    bonus = 0.0
    if in_hard:
        bonus = 0.04 * alpha + 0.03 * gamma + 0.02 * tau + 0.02 * kappa + 0.01 * mu + 0.01 * omega

    return 600 * m_norm + bonus - 0.01 * n


def _v6_sparse_long(G, props, c):
    """v6: favorise les graphes longs et peu denses (chemins, arbres)."""
    viol = violation(c, props)
    n = props.get("n", 1)
    m = props.get("m", 0)
    if viol > VIOLATION_THRESHOLD:
        return 1e6 + 500 * viol - 0.01 * n

    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)
    scale = 1 + abs(x_val) + abs(y_val) + abs(fx)
    m_norm = viol / scale

    diam = props.get("diam", 0)
    density = _density(n, m)
    # Reward large diameter, small density
    bonus = 0.05 * diam - 0.5 * density
    return 500 * m_norm + bonus - 0.01 * n


def _v7_dense_clique(G, props, c):
    """v7: favorise les graphes denses et fortement connexes."""
    viol = violation(c, props)
    n = props.get("n", 1)
    m = props.get("m", 0)
    if viol > VIOLATION_THRESHOLD:
        return 1e6 + 500 * viol - 0.01 * n

    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)
    scale = 1 + abs(x_val) + abs(y_val) + abs(fx)
    m_norm = viol / scale

    omega = props.get("omega", 0)
    kappa = props.get("kappa", 0)
    density = _density(n, m)
    # Reward large cliques, high connectivity, high density
    bonus = 0.04 * omega + 0.03 * kappa + 0.3 * density
    return 500 * m_norm + bonus - 0.015 * n


def _v8_final_mixed(G, props, c):
    """v8: combinaison optimisee structure + classe + direction."""
    viol = violation(c, props)
    n = props.get("n", 1)
    m = props.get("m", 0)
    if viol > VIOLATION_THRESHOLD:
        return 1e6 + 1e4 * viol - 0.01 * n

    x_val = get_inv(props, c["x_name"])
    y_val = get_inv(props, c["y_name"])
    fx = eval_poly(c, x_val)
    scale = 1 + abs(x_val) + abs(y_val) + abs(fx)
    m_norm = viol / scale

    diam = props.get("diam", 0)
    alpha = props.get("alpha", 0)
    gamma = props.get("gamma", 0)
    omega = props.get("omega", 0)
    mu = props.get("mu", 0)
    kappa = props.get("kappa", 0)

    # Structure
    structure = (0.025 * diam + 0.02 * alpha + 0.015 * gamma +
                 0.01 * omega + 0.008 * mu + 0.005 * kappa)

    # Class bonus
    sg = (c["subgroup"][0] if c.get("subgroup") else "connected").lower()
    class_bonus = 0.0
    if sg == "tree":
        class_bonus = 0.02 * diam
    elif sg == "claw_free":
        class_bonus = 0.02 * omega
    elif sg == "bipartite":
        class_bonus = 0.015 * mu

    # Direction
    direction = 0.04 * viol
    signal = 0.03 * abs(viol) + 0.01 * abs(x_val)

    density = _density(n, m)
    complexity = 0.018 * n + 0.003 * m + 0.06 * abs(density - 0.45)

    return 700 * m_norm + structure + class_bonus + direction + signal - complexity


SCORE_FUNCTIONS = {
    "v0_violation_pure": _v0_violation_pure,
    "v1_normalized_margin": _v1_normalized_margin,
    "v2_structure_guided": _v2_structure_guided,
    "v3_conjecture_direction": _v3_conjecture_direction,
    "v4_class_specialized": _v4_class_specialized,
    "v5_hard_invariant_guided": _v5_hard_invariant_guided,
    "v6_sparse_long": _v6_sparse_long,
    "v7_dense_clique": _v7_dense_clique,
    "v8_final_mixed": _v8_final_mixed,
}


def choose_score(G, props, c):
    """Auto-selection: evaluate all 9 functions, return best (name, value)."""
    best_name, best_val = "v0_violation_pure", -math.inf
    for name, fn in SCORE_FUNCTIONS.items():
        try:
            val = fn(G, props, c)
            if val > best_val:
                best_val = val
                best_name = name
        except Exception:
            pass
    return best_name, best_val


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
                    for w in comps[i + 1]:
                        if color.get(w, 0) != color.get(u, 0):
                            G.add_edge(u, w)
                            break

    elif sg == "planar":
        if not nx.is_planar(G):
            G = nx.minimum_spanning_tree(G)

    return G


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
            G = nx.bipartite.random_graph(half, n - half, 0.5,
                                          seed=rng.randint(0, 99999))
            if nx.is_connected(G):
                graphs.append(nx.convert_node_labels_to_integers(G))
    elif sg == "claw_free":
        for n in range(3, min(max_n + 1, 14)):
            for base in [nx.path_graph(n), nx.cycle_graph(n), nx.complete_graph(n)]:
                LG = nx.convert_node_labels_to_integers(nx.line_graph(base))
                if nx.is_connected(LG):
                    graphs.append(LG)
        for n in range(2, min(max_n + 1, 14)):
            graphs.append(nx.complete_graph(n))
            graphs.append(nx.cycle_graph(max(3, n)))
    elif sg == "planar":
        for n in range(2, min(max_n + 1, 16)):
            graphs.append(nx.path_graph(n))
            if n >= 3:
                graphs.append(nx.cycle_graph(n))
            graphs.append(nx.ladder_graph(n))
        for k in range(2, 6):
            graphs.append(nx.grid_2d_graph(k, k))
        graphs += [nx.petersen_graph(), nx.dodecahedral_graph()]
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

    for G in [nx.petersen_graph(), nx.dodecahedral_graph(),
              nx.icosahedral_graph(), nx.octahedral_graph(), nx.cubical_graph()]:
        bank.append(G)

    result = []
    for G in bank:
        G = nx.convert_node_labels_to_integers(G)
        if 1 <= G.number_of_nodes() <= max_n + 4 and nx.is_connected(G):
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

    if op == "add_leaf" and n < max_n:
        G.add_edge(rng.choice(nodes), max(nodes) + 1)
    elif op == "remove_leaf":
        leaves = [v for v in nodes if G.degree(v) == 1]
        if leaves:
            G.remove_node(rng.choice(leaves))
    elif op == "subdivide":
        edges = list(G.edges())
        if edges:
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
            new_p = rng.choice([v for v in non_leaves if v != leaf])
            old_p = list(G.neighbors(leaf))[0]
            if new_p != old_p:
                G.remove_edge(leaf, old_p)
                G.add_edge(leaf, new_p)
    return G


def mutate(G, subgroup, max_n, rng, steps=1):
    """Apply 1-3 mutation steps."""
    sg = (subgroup[0] if subgroup else "connected").lower()
    for _ in range(steps):
        G = _single_mutate(G, sg, max_n, rng)
    return G


def _single_mutate(G, sg, max_n, rng):
    G = G.copy()
    n = G.number_of_nodes()
    nodes = list(G.nodes())
    edges = list(G.edges())

    if sg == "tree":
        return _mutate_tree(G, max_n, rng)

    # Extended mutation set for Part 2
    ops = ["add_edge", "remove_edge", "toggle_edge", "add_vertex",
           "subdivide", "add_path", "add_clique", "add_twins",
           "densify_local", "sparsify_local", "attach_star",
           "attach_cycle", "double_subdivide"]
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

    elif op == "add_vertex" and n < max_n:
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

    elif op == "subdivide" and edges:
        u, v = rng.choice(edges)
        new_v = max(G.nodes()) + 1
        G.remove_edge(u, v)
        G.add_edge(u, new_v)
        G.add_edge(new_v, v)

    elif op == "double_subdivide" and edges:
        u, v = rng.choice(edges)
        nv1 = max(G.nodes()) + 1
        nv2 = nv1 + 1
        G.remove_edge(u, v)
        G.add_edge(u, nv1)
        G.add_edge(nv1, nv2)
        G.add_edge(nv2, v)

    elif op == "add_path" and n < max_n:
        attach = rng.choice(nodes)
        new_v = max(G.nodes()) + 1
        G.add_edge(attach, new_v)
        if n < max_n - 1 and rng.random() < 0.5:
            G.add_edge(new_v, new_v + 1)

    elif op == "add_clique" and n < max_n - 1:
        attach = rng.choice(nodes)
        nv1 = max(G.nodes()) + 1
        nv2 = nv1 + 1
        G.add_edge(attach, nv1)
        G.add_edge(attach, nv2)
        G.add_edge(nv1, nv2)

    elif op == "add_twins" and n < max_n:
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
        # Add edges within a random subset of 5 vertices
        subset = rng.sample(nodes, min(5, n))
        for i in range(len(subset)):
            for j in range(i + 1, len(subset)):
                if rng.random() < 0.4 and not G.has_edge(subset[i], subset[j]):
                    G.add_edge(subset[i], subset[j])

    elif op == "sparsify_local":
        # Remove 1-2 edges while preserving connectivity
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

    elif op == "attach_star" and n < max_n - 1:
        # Attach K_{1,2} to a vertex
        attach = rng.choice(nodes)
        nv1 = max(G.nodes()) + 1
        nv2 = nv1 + 1
        G.add_edge(attach, nv1)
        G.add_edge(attach, nv2)

    elif op == "attach_cycle" and n < max_n - 2:
        # Attach C_4 to a vertex
        attach = rng.choice(nodes)
        base = max(G.nodes()) + 1
        for i in range(4):
            G.add_edge(base + i, base + (i + 1) % 4)
        G.add_edge(attach, base)

    elif op == "line_graph" and n <= 10:
        LG = nx.convert_node_labels_to_integers(nx.line_graph(G))
        if nx.is_connected(LG) and 1 <= LG.number_of_nodes() <= max_n + 2:
            G = LG

    return G


def _fresh_candidate(sg, max_n, rng):
    """Generate a fresh random candidate for the class."""
    n = rng.randint(3, max(4, max_n - 2))
    if sg == "tree":
        if rng.random() < 0.5:
            seq = [rng.randint(0, n - 1) for _ in range(max(1, n - 2))]
            try:
                return nx.from_prufer_sequence(seq)
            except Exception:
                return nx.path_graph(n)
        return nx.star_graph(n)
    elif sg == "bipartite":
        a = rng.randint(1, max(2, n // 2))
        b = n - a
        G = nx.bipartite.random_graph(a, b, rng.random() * 0.7 + 0.1,
                                      seed=rng.randint(0, 99999))
        if nx.is_connected(G):
            return nx.convert_node_labels_to_integers(G)
        return nx.complete_bipartite_graph(a, b)
    elif sg == "claw_free":
        # Line graph of random graph
        base = nx.gnm_random_graph(n, rng.randint(n, 2 * n), seed=rng.randint(0, 99999))
        LG = nx.convert_node_labels_to_integers(nx.line_graph(base))
        if nx.is_connected(LG):
            return LG
        return nx.complete_graph(n)
    else:
        p = rng.random() * 0.5 + 0.2
        G = nx.gnp_random_graph(n, p, seed=rng.randint(0, 99999))
        if nx.is_connected(G):
            return nx.convert_node_labels_to_integers(G)
        return nx.path_graph(n)


# ══════════════════════════════════════════════════════ MAIN SOLVER ═══════════

def _make_result(c, G, props, viol, cost, status, score_variant, exec_time=None):
    g6 = ""
    order = size = 0
    x_val = y_val = 0.0
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


def solve_one(c, initial_bank, seconds, max_n, rng):
    t0 = time.time()
    subgroup = c["subgroup"]
    sg = (subgroup[0] if subgroup else "connected").lower()

    class_bank = [G for G in initial_bank if is_in_class(G, subgroup)]
    if len(class_bank) < 5:
        class_bank = _class_specific(sg, max_n, rng)
    class_bank = class_bank[:300]

    best_score = -math.inf
    best_name = "v8_final_mixed"
    best_G = None
    best_props = None
    best_viol = -math.inf

    population = []

    # Evaluate initial bank
    for G in class_bank:
        if time.time() - t0 > seconds * 0.25:
            break
        try:
            props = exact_parameters(G)
            viol = violation(c, props)
            if viol > VIOLATION_THRESHOLD:
                return _make_result(c, G, props, viol, time.time() - t0, "FOUND", "v8_final_mixed")
            name, s = choose_score(G, props, c)
            population.append((s, G.copy(), props, name))
            if s > best_score:
                best_score, best_G, best_props, best_viol, best_name = s, G.copy(), props, viol, name
        except Exception:
            pass

    population.sort(key=lambda x: x[0], reverse=True)
    population = population[:140]

    # Population selection weights for top-k brackets
    TOP_BRACKETS = [12, 20, 35, 60]

    while time.time() - t0 < seconds:
        r = rng.random()

        # Select parent
        if r < 0.78 and population:
            k = min(len(population), rng.choice(TOP_BRACKETS))
            idx = int(rng.random() ** 1.5 * k)
            parent = population[idx][1]
        elif r < 0.93:
            try:
                parent = _fresh_candidate(sg, max_n, rng)
                parent = nx.convert_node_labels_to_integers(parent)
            except Exception:
                parent = rng.choice(class_bank) if class_bank else nx.path_graph(4)
        else:
            parent = rng.choice(class_bank) if class_bank else nx.path_graph(4)

        # Number of mutation steps: 1 (70%), 2 (20%), 3 (10%)
        r2 = rng.random()
        steps = 1 if r2 < 0.70 else (2 if r2 < 0.90 else 3)

        try:
            H = mutate(parent, subgroup, max_n, rng, steps)
            H = repair_graph(H, subgroup, rng)
            if H.number_of_nodes() < 1 or H.number_of_nodes() > max_n + 3:
                continue

            props = exact_parameters(H)
            viol = violation(c, props)

            if viol > VIOLATION_THRESHOLD:
                _, s = choose_score(H, props, c)
                name = "v8_final_mixed"
                for fn_name, fn in SCORE_FUNCTIONS.items():
                    try:
                        if fn(H, props, c) == s:
                            name = fn_name
                            break
                    except Exception:
                        pass
                return _make_result(c, H, props, viol, time.time() - t0, "FOUND", name)

            fn_name, s = choose_score(H, props, c)
            if len(population) < 140 or s > population[-1][0]:
                population.append((s, H.copy(), props, fn_name))
                population.sort(key=lambda x: x[0], reverse=True)
                population = population[:140]

            if s > best_score:
                best_score, best_G, best_props, best_viol, best_name = s, H.copy(), props, viol, fn_name
        except Exception:
            pass

    exec_time = time.time() - t0
    return _make_result(c, best_G, best_props,
                        best_viol if best_viol != -math.inf else -999,
                        120, "NOT_FOUND", best_name, exec_time)


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
        "conjecture_id", "status", "validation_status", "score_variant",
        "counterexample_g6", "order", "size", "counter_x", "counter_y",
        "violation_margin", "cost", "execution_time_seconds",
    ]
    total_cost = 0
    found = 0

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, c in enumerate(conjectures):
            label = f"[{i+1}/{len(conjectures)}] #{c['id']}"
            print(f"{label:20s} {c['text'][:55]}", end=" ... ", flush=True)
            try:
                res = solve_one(c, bank, seconds, max_order, random.Random(seed + c["id"]))
            except Exception as e:
                print(f"ERROR: {e}")
                res = _make_result(c, None, None, -999, 120, "NOT_FOUND", "v8_final_mixed")
            total_cost += res["cost"]
            if res["status"] == "FOUND":
                found += 1
                print(f"FOUND  t={res['cost']:.2f}s  variant={res['score_variant']}  "
                      f"margin={res['violation_margin']}")
            else:
                print("NOT FOUND  (cost=120)")
            writer.writerow(res)
            f.flush()

    print(f"\n{'='*60}")
    print(f"Conjectures refutees : {found} / {len(conjectures)}")
    print(f"Score total          : {total_cost:.2f} s")
    print(f"Resultats            : {output_path}")
    print(f"{'='*60}")


def main():
    p = argparse.ArgumentParser(
        description="GraphBench Partie 2 - FunSearch-inspired (portfolio de 9 scores)")
    p.add_argument("--input",     required=True,         help="benchmark.xlsx")
    p.add_argument("--output",    required=True,         help="fichier CSV resultats")
    p.add_argument("--seconds",   type=float, default=60.0,  help="limite temps / conjecture")
    p.add_argument("--max-order", type=int,   default=18,    help="ordre max des graphes")
    p.add_argument("--seed",      type=int,   default=42,    help="graine aleatoire")
    args = p.parse_args()
    solve(args.input, args.output, args.seconds, args.max_order, args.seed)


if __name__ == "__main__":
    main()
