"""Graph-theoretic invariants and bipartite entanglement of graph states.

For a graph state |G> and a bipartition A|B of its qubits, let Gamma_AB be
the |A| x |B| adjacency matrix between A and B, and let

    r(A) = rank of Gamma_AB over GF(2)

be the cut-rank of A. The reduced state rho_A has exactly 2**r(A) non-zero
eigenvalues, all equal to 2**(-r(A)) [Hein, Eisert & Briegel, Phys. Rev. A
69, 062311 (2004)]. The Schmidt rank across A|B is therefore 2**r(A), and the
entanglement entropy across A|B is r(A) ebits.

Schmidt ranks do not change under local unitaries. The distribution of
Schmidt ranks over all bipartitions is therefore the same for every graph in
an LC orbit, including after relabelling the qubits. Its maximum gives a lower
bound on the Schmidt measure: E_S(|G>) >= log2(max_A 2**r(A)).

All functions take a simple undirected ``networkx.Graph``. None of them
depends on Tk, so they can be used from scripts and worker processes.
"""

from collections import Counter
import warnings

import networkx as nx
import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
#  GF(2) linear algebra
# ─────────────────────────────────────────────────────────────────────────────

def gf2_rank(M) -> int:
    """Rank of a binary matrix over GF(2), by Gaussian elimination."""
    A = np.asarray(M)
    if A.size == 0:
        return 0
    A = (A.astype(np.int64) & 1).astype(np.uint8)
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if A[i, c]), None)
        if pivot is None:
            continue
        if pivot != r:
            A[[r, pivot]] = A[[pivot, r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i, :] ^= A[r, :]
        r += 1
        if r == rows:
            break
    return r


def _rank_of_bitrows(rows) -> int:
    """GF(2) rank of row vectors stored as Python integers (bit i = column i)."""
    pivots = {}
    for row in rows:
        while row:
            top = row.bit_length() - 1
            if top in pivots:
                row ^= pivots[top]
            else:
                pivots[top] = row
                break
    return len(pivots)


def _set_bits(x):
    """Indices of the set bits of ``x``, lowest first."""
    while x:
        low = x & -x
        yield low.bit_length() - 1
        x ^= low


def _neighbour_masks(G, nodes):
    """Adjacency of ``G`` as one bit mask per vertex, in the order of ``nodes``."""
    index = {v: i for i, v in enumerate(nodes)}
    masks = [0] * len(nodes)
    for u, v in G.edges():
        if u == v:
            raise ValueError("graph states are defined on graphs without self-loops")
        masks[index[u]] |= 1 << index[v]
        masks[index[v]] |= 1 << index[u]
    return masks


# ─────────────────────────────────────────────────────────────────────────────
#  Bipartitions and Schmidt ranks
# ─────────────────────────────────────────────────────────────────────────────

def cut_matrix(G: nx.Graph, A_nodes, B_nodes) -> np.ndarray:
    """Adjacency matrix Gamma_AB between the two sides of a bipartition.

    Entry (i, j) is 1 if and only if ``A_nodes[i]`` and ``B_nodes[j]`` are
    adjacent. Edges inside A or inside B do not appear.

    Raises
    ------
    ValueError
        If A and B overlap, contain duplicates, or do not cover all nodes.
    """
    A_nodes = list(A_nodes)
    B_nodes = list(B_nodes)
    a, b = set(A_nodes), set(B_nodes)
    if len(a) != len(A_nodes) or len(b) != len(B_nodes):
        raise ValueError("a side of the bipartition contains a node twice")
    if a & b:
        raise ValueError("the two sides of the bipartition overlap")
    if a | b != set(G.nodes):
        raise ValueError("the two sides of the bipartition must cover every node")

    idx_a = {u: i for i, u in enumerate(A_nodes)}
    idx_b = {v: j for j, v in enumerate(B_nodes)}
    M = np.zeros((len(A_nodes), len(B_nodes)), dtype=np.uint8)
    for u, v in G.edges():
        if u in idx_a and v in idx_b:
            M[idx_a[u], idx_b[v]] = 1
        elif v in idx_a and u in idx_b:
            M[idx_a[v], idx_b[u]] = 1
    return M


def cut_rank(G: nx.Graph, A_nodes) -> int:
    """Cut-rank r(A): GF(2) rank of Gamma_AB with B the complement of A.

    Equals the entanglement entropy of |G> across A|B in ebits.
    """
    A = list(A_nodes)
    a = set(A)
    if not a <= set(G.nodes):
        raise ValueError("A contains nodes that are not in G")
    B = [v for v in G.nodes if v not in a]
    return gf2_rank(cut_matrix(G, A, B))


def schmidt_rank(G: nx.Graph, A_nodes) -> int:
    """Schmidt rank 2**r(A) of the graph state |G> across A|B."""
    return 2 ** cut_rank(G, A_nodes)


def schmidt_rank_profile(G: nx.Graph, max_nodes: int | None = None):
    """Schmidt ranks of |G> over all bipartitions.

    Every unordered bipartition {A, B} with A and B non-empty is counted once,
    so the counts add up to 2**(n-1) - 1.

    Parameters
    ----------
    G : networkx.Graph
        Simple undirected graph on n vertices.
    max_nodes : int or None
        Refuse graphs with more than ``max_nodes`` vertices. The cost grows as
        2**(n-1).

    Returns
    -------
    max_rank : int
        Largest Schmidt rank over all bipartitions (1 when n < 2).
    distribution : dict
        Maps each Schmidt rank 2**r to the number of bipartitions with that
        rank, in increasing order of rank.
    """
    nodes = list(G.nodes)
    n = len(nodes)
    if max_nodes is not None and n > max_nodes:
        raise ValueError(
            f"{n} vertices exceed the limit of {max_nodes} for exhaustive enumeration"
        )
    if n < 2:
        return 1, {}

    masks = _neighbour_masks(G, nodes)
    full = (1 << n) - 1
    counts = Counter()
    # Fixing vertex 0 on side A lists each unordered bipartition exactly once.
    for m in range(1 << (n - 1)):
        A = (m << 1) | 1
        if A == full:
            continue
        side = A if A.bit_count() <= n // 2 else full ^ A
        other = full ^ side
        counts[_rank_of_bitrows(masks[i] & other for i in _set_bits(side))] += 1

    distribution = {2 ** r: counts[r] for r in sorted(counts)}
    return max(distribution), distribution


# ─────────────────────────────────────────────────────────────────────────────
#  Colouring
# ─────────────────────────────────────────────────────────────────────────────

def greedy_coloring(G: nx.Graph) -> dict:
    """Greedy colouring with the saturation-largest-first (DSATUR) strategy.

    Returns a mapping node -> colour index. The number of colours used is an
    upper bound on the chromatic number.
    """
    return nx.coloring.greedy_color(G, strategy="saturation_largest_first")


def is_k_colorable(G: nx.Graph, k: int, max_steps: int = 200_000):
    """Decide whether ``G`` admits a proper colouring with ``k`` colours.

    Exact backtracking search with DSATUR vertex ordering.

    Returns
    -------
    status : bool or None
        True or False when the search finished, None when it stopped after
        ``max_steps`` steps without a decision.
    coloring : dict
        A proper colouring node -> colour index when ``status`` is True,
        otherwise an empty dict.
    """
    nodes = list(G.nodes)
    if not nodes:
        return True, {}
    if k < 1:
        return False, {}
    if k == 1:
        if G.number_of_edges() == 0:
            return True, {v: 0 for v in nodes}
        return False, {}

    colors = {}
    steps = 0

    def saturation(v):
        return len({colors[u] for u in G.neighbors(v) if u in colors})

    def next_vertex():
        uncolored = [v for v in nodes if v not in colors]
        return max(uncolored, key=lambda v: (saturation(v), G.degree[v]))

    def search():
        nonlocal steps
        if len(colors) == len(nodes):
            return True
        steps += 1
        if steps > max_steps:
            return None
        v = next_vertex()
        used = {colors[u] for u in G.neighbors(v) if u in colors}
        for c in range(k):
            if c in used:
                continue
            colors[v] = c
            found = search()
            if found is not False:
                return found
            del colors[v]
        return False

    status = search()
    return status, (dict(colors) if status else {})


# ─────────────────────────────────────────────────────────────────────────────
#  Structural statistics
# ─────────────────────────────────────────────────────────────────────────────

def _largest_component(G: nx.Graph) -> nx.Graph:
    """Largest connected component of ``G`` (``G`` itself if connected)."""
    if G.number_of_nodes() == 0 or nx.is_connected(G):
        return G
    return G.subgraph(max(nx.connected_components(G), key=len)).copy()


def _assortativity(G: nx.Graph):
    """Degree assortativity, or None where it is undefined (e.g. regular graphs)."""
    if G.number_of_edges() == 0:
        return None
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        value = nx.degree_assortativity_coefficient(G)
    return None if np.isnan(value) else float(value)


def graph_stats(G: nx.Graph, schmidt_max_nodes: int = 16) -> dict:
    """Structural invariants of ``G`` and bipartite entanglement of |G>.

    Parameters
    ----------
    G : networkx.Graph
        Simple undirected graph.
    schmidt_max_nodes : int
        Schmidt ranks are computed only up to this many vertices; above it the
        Schmidt entries are None.

    Returns
    -------
    dict
        ``nodes``, ``edges``, ``density``, ``connected``, ``components``,
        ``diameter_giant`` and ``avg_path_len_giant`` (largest component),
        ``clustering``, ``assortativity``, ``bipartite``, ``bipartition_sizes``
        (connected bipartite graphs only), ``planar``, ``clique_number``,
        ``greedy_chromatic``, ``chromatic_lower_bound``,
        ``chromatic_upper_bound``, ``deg_min``, ``deg_avg``, ``deg_max``,
        ``triangles``, ``schmidt_rank_max`` and ``schmidt_distribution``
        (see :func:`schmidt_rank_profile`). Entries that are undefined for the
        given graph are None.
    """
    n = G.number_of_nodes()
    m = G.number_of_edges()
    degrees = [d for _, d in G.degree()]

    giant = _largest_component(G)
    connected = nx.is_connected(G) if n > 0 else None
    bipartite = nx.is_bipartite(G)

    bipartition_sizes = None
    if bipartite and connected:
        side_a, side_b = nx.bipartite.sets(G)
        bipartition_sizes = (len(side_a), len(side_b))

    clique_number = max((len(c) for c in nx.find_cliques(G)), default=0)
    greedy = greedy_coloring(G)
    greedy_chi = len(set(greedy.values()))

    schmidt_max = None
    schmidt_dist = None
    if 0 < n <= schmidt_max_nodes:
        schmidt_max, schmidt_dist = schmidt_rank_profile(G)

    return {
        "nodes": n,
        "edges": m,
        "density": float(nx.density(G)) if n > 1 else 0.0,
        "connected": connected,
        "components": nx.number_connected_components(G) if n > 0 else 0,
        "diameter_giant": nx.diameter(giant) if n > 0 else None,
        "avg_path_len_giant": (
            float(nx.average_shortest_path_length(giant)) if n > 0 else None
        ),
        "clustering": float(nx.average_clustering(G)) if n > 0 else 0.0,
        "assortativity": _assortativity(G),
        "bipartite": bipartite,
        "bipartition_sizes": bipartition_sizes,
        "planar": nx.check_planarity(G)[0],
        "clique_number": clique_number,
        "greedy_chromatic": greedy_chi,
        "chromatic_lower_bound": clique_number,
        "chromatic_upper_bound": greedy_chi,
        "deg_min": min(degrees, default=0),
        "deg_avg": round(float(np.mean(degrees)), 3) if degrees else 0.0,
        "deg_max": max(degrees, default=0),
        "triangles": sum(nx.triangles(G).values()) // 3,
        "schmidt_rank_max": schmidt_max,
        "schmidt_distribution": schmidt_dist,
    }


ORBIT_TABLE_COLUMNS = (
    "graph6", "n", "m", "connected", "diam(Gc)", "planar", "bipartite",
    "omega", "chi_greedy", "Schmidt_max",
)


def orbit_table_row(g6: str) -> tuple:
    """One row of the orbit statistics table, in the order of ORBIT_TABLE_COLUMNS.

    Module-level so that it can be sent to worker processes.
    """
    s = graph_stats(nx.from_graph6_bytes(g6.encode("ascii")))
    return (
        g6, s["nodes"], s["edges"], s["connected"], s["diameter_giant"],
        s["planar"], s["bipartite"], s["clique_number"], s["greedy_chromatic"],
        s["schmidt_rank_max"],
    )
