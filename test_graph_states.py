"""Tests for Graph State Orbit Explorer.

Run from the project root with::

    python -m unittest -v

The tests check the code against independent references:

* canonical labelling against the networkx graph atlas (all graphs with up to
  seven vertices, pairwise non-isomorphic);
* local complementation against a direct edge-toggling implementation;
* LC-orbit counts against the classification of Hein, Eisert & Briegel,
  Phys. Rev. A 69, 062311 (2004): 1, 1, 2, 4, 11, 26 orbits of connected
  graphs on 2 to 7 vertices;
* Schmidt ranks against the rank of the reshaped state vector of |G>;
* exact colourability against brute force.

The GUI tests need a display. They are skipped when none is available; on a
headless machine run them with ``xvfb-run python -m unittest -v``.
"""

from collections import deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from itertools import combinations, product
import multiprocessing
import random
import time
import unittest
import warnings
from unittest import mock

import networkx as nx
import numpy as np

from src import graphOperations
from src.graphAnalytics import (
    ORBIT_TABLE_COLUMNS,
    cut_matrix,
    cut_rank,
    gf2_rank,
    graph_stats,
    greedy_coloring,
    is_k_colorable,
    orbit_table_row,
    schmidt_rank,
    schmidt_rank_profile,
)
from src.toolset import canonical_relabel_n


# ─────────────────────────────────────────────────────────────────────────────
#  Reference implementations and helpers
# ─────────────────────────────────────────────────────────────────────────────

def canonical_g6(graph):
    return canonical_relabel_n(graph).decode()


def graph_from_g6(g6):
    if isinstance(g6, bytes):
        g6 = g6.decode()
    return nx.from_graph6_bytes(g6.replace(">>graph6<<", "").strip().encode())


def reference_local_complementation(graph, vertex):
    """LC by toggling every pair of neighbours of ``vertex``."""
    result = graph.copy()
    for u, v in combinations(list(graph.neighbors(vertex)), 2):
        if result.has_edge(u, v):
            result.remove_edge(u, v)
        else:
            result.add_edge(u, v)
    return result


def reference_lc_operations(graph):
    seen = set()
    operations = []
    for vertex in graph.nodes:
        if graph.degree[vertex] < 2:
            continue
        target = canonical_g6(reference_local_complementation(graph, vertex))
        if target not in seen:
            seen.add(target)
            operations.append((target, vertex))
    return operations


def reference_orbit(graph):
    start = canonical_g6(graph)
    queue = deque([start])
    seen = set()
    orbit = []
    neighbors = {}
    while queue:
        current = queue.popleft()
        if current in seen:
            continue
        seen.add(current)
        orbit.append(current)
        operations = reference_lc_operations(graph_from_g6(current))
        neighbors[current] = operations
        for target, _vertex in operations:
            if target not in seen:
                queue.append(target)
    return orbit, neighbors


def state_vector(graph):
    """|G> = prod_{(i,j) in E} CZ_ij |+>^n as an n-index tensor."""
    nodes = list(graph.nodes)
    n = len(nodes)
    index = {v: i for i, v in enumerate(nodes)}
    bits = (np.arange(2 ** n)[:, None] >> np.arange(n)[::-1]) & 1
    phase = np.zeros(2 ** n, dtype=int)
    for u, v in graph.edges():
        phase += bits[:, index[u]] * bits[:, index[v]]
    return ((-1.0) ** phase / 2 ** (n / 2)).reshape([2] * n), index


def numerical_schmidt_rank(psi, index, side):
    a = [index[v] for v in side]
    b = [i for i in range(psi.ndim) if i not in a]
    matrix = np.transpose(psi, a + b).reshape(2 ** len(a), -1)
    return int(np.linalg.matrix_rank(matrix, tol=1e-9))


def brute_force_colorable(graph, k):
    nodes = list(graph.nodes)
    position = {v: i for i, v in enumerate(nodes)}
    return any(
        all(colors[position[u]] != colors[position[v]] for u, v in graph.edges())
        for colors in product(range(k), repeat=len(nodes))
    )


class DummyProgress:
    """Stand-in for tqdm_tk that records what the progress bar was told."""

    def __init__(self, *args, **kwargs):
        self.total = kwargs.get("total")
        self.n = 0
        self.descriptions = []
        self.closed = False
        self.refreshed = False

    def set_description(self, text):
        self.descriptions.append(text)

    def update(self, n=1):
        self.n += n

    def refresh(self):
        self.refreshed = True

    def close(self):
        self.closed = True


@contextmanager
def headless_progress():
    bars = []

    def make(*args, **kwargs):
        bar = DummyProgress(*args, **kwargs)
        bars.append(bar)
        return bar

    with mock.patch.object(graphOperations, "tqdm_tk", side_effect=make):
        yield bars


def explore(seeds):
    with headless_progress():
        return graphOperations.exploreOrbit(seeds, master=None)


SMALL_GRAPHS = {
    "empty_4": nx.empty_graph(4),
    "path_4": nx.path_graph(4),
    "cycle_4": nx.cycle_graph(4),
    "star_4": nx.star_graph(3),
    "complete_4": nx.complete_graph(4),
    "paw": nx.Graph([(0, 1), (1, 2), (2, 0), (2, 3)]),
    "house": nx.house_graph(),
    "path_3_plus_isolate": nx.disjoint_union(nx.path_graph(3), nx.empty_graph(1)),
}


# ─────────────────────────────────────────────────────────────────────────────
#  Canonical labelling
# ─────────────────────────────────────────────────────────────────────────────

class CanonicalLabellingTests(unittest.TestCase):
    def test_distinguishes_every_graph_with_up_to_seven_vertices(self):
        atlas = nx.graph_atlas_g()
        strings = {canonical_g6(g) for g in atlas}
        self.assertEqual(len(strings), len(atlas))

    def test_invariant_under_relabelling(self):
        rng = random.Random(7)
        for graph in rng.sample(nx.graph_atlas_g()[1:], 100):
            reference = canonical_g6(graph)
            nodes = list(graph.nodes)
            shuffled = nodes[:]
            rng.shuffle(shuffled)
            variants = [
                nx.relabel_nodes(graph, dict(zip(nodes, shuffled))),
                nx.relabel_nodes(graph, {v: 10 * v + 3 for v in nodes}),
                nx.relabel_nodes(graph, {v: f"q{v}" for v in nodes}),
            ]
            for variant in variants:
                self.assertEqual(canonical_g6(variant), reference)

    def test_output_format_and_empty_graph(self):
        self.assertEqual(canonical_relabel_n(nx.Graph()), b"?")
        self.assertEqual(graph_from_g6(canonical_relabel_n(nx.Graph())).number_of_nodes(), 0)
        g6 = canonical_relabel_n(nx.petersen_graph())
        self.assertNotIn(b"\n", g6)
        self.assertFalse(g6.startswith(b">>"))

    def test_rejects_graphs_that_are_not_simple(self):
        with self.assertRaises(ValueError):
            canonical_relabel_n(nx.Graph([(0, 0), (0, 1)]))
        with self.assertRaises(TypeError):
            canonical_relabel_n(nx.DiGraph([(0, 1)]))
        with self.assertRaises(TypeError):
            canonical_relabel_n(nx.MultiGraph([(0, 1)]))


# ─────────────────────────────────────────────────────────────────────────────
#  Local complementation
# ─────────────────────────────────────────────────────────────────────────────

class LocalComplementationTests(unittest.TestCase):
    def test_centre_of_path_gives_triangle_and_leaf_does_nothing(self):
        path = nx.path_graph(3)
        centre = graph_from_g6(graphOperations.local_complementation_at_vertex(path, 1))
        leaf = graph_from_g6(graphOperations.local_complementation_at_vertex(path, 0))
        self.assertTrue(nx.is_isomorphic(centre, nx.complete_graph(3)))
        self.assertEqual(sorted(leaf.edges), sorted(path.edges))

    def test_matches_toggle_rule_at_every_vertex(self):
        rng = random.Random(11)
        graphs = list(SMALL_GRAPHS.values()) + [
            nx.gnp_random_graph(rng.randint(2, 8), 0.5, seed=rng.randrange(10 ** 6))
            for _ in range(40)
        ]
        for graph in graphs:
            for vertex in graph.nodes:
                result = graph_from_g6(graphOperations.local_complementation_at_vertex(graph, vertex))
                expected = reference_local_complementation(graph, vertex)
                self.assertEqual(sorted(result.edges), sorted(expected.edges))

    def test_is_an_involution(self):
        graph = nx.petersen_graph()
        for vertex in graph.nodes:
            once = graph_from_g6(graphOperations.local_complementation_at_vertex(graph, vertex))
            twice = graph_from_g6(graphOperations.local_complementation_at_vertex(once, vertex))
            self.assertEqual(sorted(twice.edges), sorted(graph.edges))

    def test_vertex_position_is_validated(self):
        path = nx.path_graph(3)
        for bad in (-1, 3):
            with self.assertRaises(IndexError):
                graphOperations.local_complementation_at_vertex(path, bad)
        with self.assertRaises(TypeError):
            graphOperations.local_complementation_at_vertex(path, 1.0)

    def test_edge_weights_are_ignored(self):
        triangle = nx.Graph()
        triangle.add_edge(0, 1)
        triangle.add_edge(0, 2)
        triangle.add_edge(1, 2, weight=2)
        result = graph_from_g6(graphOperations.local_complementation_at_vertex(triangle, 0))
        self.assertEqual(sorted(result.edges), [(0, 1), (0, 2)])

    def test_neighbour_enumeration_matches_reference(self):
        for name, graph in SMALL_GRAPHS.items():
            with self.subTest(graph=name):
                expected = reference_lc_operations(graph)
                targets, operations = graphOperations.localComplementationOperationUI(graph)
                self.assertEqual(operations, expected)
                self.assertEqual(targets, [t for t, _v in expected])
                self.assertEqual(graphOperations.localComplementationOperation(graph), targets)


# ─────────────────────────────────────────────────────────────────────────────
#  Orbit exploration
# ─────────────────────────────────────────────────────────────────────────────

def connected_orbits(max_nodes):
    """All LC orbits of connected graphs on 2..max_nodes vertices, by n."""
    orbits = {}
    for n in range(2, max_nodes + 1):
        seen = set()
        orbits[n] = []
        for graph in nx.graph_atlas_g():
            if graph.number_of_nodes() != n or not nx.is_connected(graph):
                continue
            g6 = canonical_g6(graph)
            if g6 in seen:
                continue
            orbit, _neighbors = explore([g6])
            seen.update(orbit)
            orbits[n].append(orbit)
    return orbits


class OrbitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.orbits = connected_orbits(7)

    def test_number_of_orbits_matches_literature(self):
        counts = [len(self.orbits[n]) for n in range(2, 8)]
        self.assertEqual(counts, [1, 1, 2, 4, 11, 26])

    def test_orbits_partition_the_connected_graphs(self):
        for n, orbits in self.orbits.items():
            members = [g6 for orbit in orbits for g6 in orbit]
            connected = {canonical_g6(g) for g in nx.graph_atlas_g()
                         if g.number_of_nodes() == n and nx.is_connected(g)}
            self.assertEqual(len(members), len(set(members)))
            self.assertEqual(set(members), connected)

    def test_schmidt_distribution_is_constant_on_every_orbit(self):
        for orbits in self.orbits.values():
            for orbit in orbits:
                profiles = {repr(schmidt_rank_profile(graph_from_g6(g6))) for g6 in orbit}
                self.assertEqual(len(profiles), 1, orbit)

    def test_complete_graph_and_star_form_one_orbit(self):
        for n in range(3, 8):
            orbit, _ = explore([nx.to_graph6_bytes(nx.complete_graph(n), header=False)])
            self.assertEqual(set(orbit), {canonical_g6(nx.complete_graph(n)),
                                          canonical_g6(nx.star_graph(n - 1))})

    def test_matches_reference_search(self):
        for name, graph in SMALL_GRAPHS.items():
            with self.subTest(graph=name):
                expected_orbit, expected_neighbors = reference_orbit(graph)
                with headless_progress() as bars:
                    orbit, neighbors = graphOperations.exploreOrbit(
                        [nx.to_graph6_bytes(graph, header=False)], master=None)
                self.assertEqual(orbit, expected_orbit)
                self.assertEqual(neighbors, expected_neighbors)
                self.assertEqual(bars[0].n, bars[0].total)
                self.assertTrue(bars[0].refreshed and bars[0].closed)

    def test_seed_formats_and_isomorphic_seeds(self):
        path = nx.path_graph(4)
        relabelled = nx.relabel_nodes(path, {0: 3, 1: 1, 2: 0, 3: 2})
        expected, _ = reference_orbit(path)
        for seed in (nx.to_graph6_bytes(relabelled, header=False),
                     nx.to_graph6_bytes(relabelled).decode(),
                     canonical_g6(path)):
            orbit, _ = explore([seed])
            self.assertEqual(orbit, expected)

    def test_runs_without_a_progress_window(self):
        expected, _ = reference_orbit(nx.cycle_graph(5))
        with mock.patch.object(graphOperations, "tqdm_tk") as window:
            orbit, _ = graphOperations.exploreOrbit(
                [canonical_g6(nx.cycle_graph(5))], progress=False)
        window.assert_not_called()
        self.assertEqual(orbit, expected)

    def test_input_is_not_modified(self):
        seeds = [nx.to_graph6_bytes(nx.path_graph(4), header=False)]
        copy = list(seeds)
        explore(seeds)
        self.assertEqual(seeds, copy)

    def test_seeds_from_one_orbit_and_from_two_orbits(self):
        path = nx.to_graph6_bytes(nx.path_graph(3), header=False)
        triangle = nx.to_graph6_bytes(nx.complete_graph(3), header=False)
        orbit, _ = explore([path, triangle])
        self.assertEqual(set(orbit), {canonical_g6(nx.path_graph(3)),
                                      canonical_g6(nx.complete_graph(3))})

        k5 = nx.to_graph6_bytes(nx.complete_graph(5), header=False)
        p4 = nx.to_graph6_bytes(nx.path_graph(4), header=False)
        orbit, neighbors = explore([p4, k5])
        expected = set(reference_orbit(nx.path_graph(4))[0]) | set(reference_orbit(nx.complete_graph(5))[0])
        self.assertEqual(set(orbit), expected)
        self.assertEqual(set(neighbors), expected)


# ─────────────────────────────────────────────────────────────────────────────
#  Schmidt ranks and graph statistics
# ─────────────────────────────────────────────────────────────────────────────

class SchmidtRankTests(unittest.TestCase):
    def test_gf2_rank(self):
        dependent = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=np.uint8)
        self.assertEqual(gf2_rank(dependent), 2)
        self.assertEqual(gf2_rank(np.eye(4, dtype=np.uint8)), 4)
        self.assertEqual(gf2_rank(np.zeros((3, 2), dtype=np.uint8)), 0)
        self.assertEqual(gf2_rank(np.zeros((0, 3), dtype=np.uint8)), 0)

    def test_cut_matrix_and_ranks_of_path(self):
        path = nx.path_graph(4)
        np.testing.assert_array_equal(cut_matrix(path, [0, 2], [1, 3]),
                                      np.array([[1, 0], [1, 1]], dtype=np.uint8))
        self.assertEqual(cut_rank(path, [0, 2]), 2)
        self.assertEqual(schmidt_rank(path, [0, 2]), 4)
        self.assertEqual(schmidt_rank(path, [0, 1]), 2)

    def test_cut_matrix_rejects_invalid_bipartitions(self):
        path = nx.path_graph(4)
        with self.assertRaises(ValueError):
            cut_matrix(path, [0, 1], [1, 2, 3])
        with self.assertRaises(ValueError):
            cut_matrix(path, [0], [1, 2])
        with self.assertRaises(ValueError):
            cut_matrix(path, [0, 0], [1, 2, 3])
        with self.assertRaises(ValueError):
            cut_rank(path, [7])

    def test_schmidt_rank_matches_state_vector(self):
        rng = random.Random(5)
        for _ in range(25):
            n = rng.randint(2, 7)
            graph = nx.gnp_random_graph(n, rng.choice([0.3, 0.6, 0.9]), seed=rng.randrange(10 ** 6))
            psi, index = state_vector(graph)
            for k in range(1, n // 2 + 1):
                for side in combinations(graph.nodes, k):
                    self.assertEqual(schmidt_rank(graph, side),
                                     numerical_schmidt_rank(psi, index, side))

    def test_profiles_of_standard_states(self):
        self.assertEqual(schmidt_rank_profile(nx.path_graph(4)), (4, {2: 5, 4: 2}))
        self.assertEqual(schmidt_rank_profile(nx.complete_graph(3)), (2, {2: 3}))
        self.assertEqual(schmidt_rank_profile(nx.complete_graph(6)), (2, {2: 31}))
        self.assertEqual(schmidt_rank_profile(nx.cycle_graph(6)), (8, {2: 6, 4: 19, 8: 6}))
        self.assertEqual(schmidt_rank_profile(nx.empty_graph(1)), (1, {}))

    def test_maximum_is_not_taken_across_the_colour_classes(self):
        k33 = nx.complete_bipartite_graph(3, 3)
        self.assertEqual(schmidt_rank(k33, [0, 1, 2]), 2)
        self.assertEqual(schmidt_rank_profile(k33)[0], 4)

    def test_profile_counts_each_bipartition_once(self):
        rng = random.Random(9)
        for n in range(2, 10):
            graph = nx.gnp_random_graph(n, 0.5, seed=rng.randrange(10 ** 6))
            max_rank, dist = schmidt_rank_profile(graph)
            self.assertEqual(sum(dist.values()), 2 ** (n - 1) - 1)
            self.assertEqual(max_rank, max(dist))
            nodes = list(graph.nodes)
            brute = {}
            for k in range(1, n):
                for side in combinations(nodes, k):
                    if nodes[0] in side:
                        rank = schmidt_rank(graph, side)
                        brute[rank] = brute.get(rank, 0) + 1
            self.assertEqual(dist, dict(sorted(brute.items())))

    def test_size_limit(self):
        with self.assertRaises(ValueError):
            schmidt_rank_profile(nx.cycle_graph(10), max_nodes=9)


class GraphStatsTests(unittest.TestCase):
    def test_triangle(self):
        stats = graph_stats(nx.complete_graph(3))
        self.assertEqual(stats["nodes"], 3)
        self.assertEqual(stats["edges"], 3)
        self.assertEqual(stats["clique_number"], 3)
        self.assertEqual(stats["greedy_chromatic"], 3)
        self.assertEqual(stats["triangles"], 1)
        self.assertEqual(stats["schmidt_rank_max"], 2)
        self.assertEqual(stats["schmidt_distribution"], {2: 3})
        self.assertFalse(stats["bipartite"])
        self.assertTrue(stats["planar"])

    def test_agrees_with_networkx(self):
        rng = random.Random(3)
        for _ in range(60):
            graph = nx.gnp_random_graph(rng.randint(1, 9), rng.random(), seed=rng.randrange(10 ** 6))
            stats = graph_stats(graph)
            giant = graph.subgraph(max(nx.connected_components(graph), key=len))
            self.assertEqual(stats["edges"], graph.number_of_edges())
            self.assertEqual(stats["components"], nx.number_connected_components(graph))
            self.assertEqual(stats["diameter_giant"], nx.diameter(giant))
            self.assertEqual(stats["clique_number"], max(len(c) for c in nx.find_cliques(graph)))
            self.assertEqual(stats["triangles"], sum(nx.triangles(graph).values()) // 3)
            self.assertEqual(stats["planar"], nx.check_planarity(graph)[0])
            self.assertEqual(stats["bipartite"], nx.is_bipartite(graph))

    def test_edge_cases_without_warnings(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            null = graph_stats(nx.Graph())
            regular = graph_stats(nx.cycle_graph(5))
            split = graph_stats(nx.disjoint_union(nx.path_graph(3), nx.path_graph(2)))
        self.assertIsNone(null["connected"])
        self.assertIsNone(null["schmidt_rank_max"])
        self.assertIsNone(regular["assortativity"])
        self.assertTrue(split["bipartite"])
        self.assertIsNone(split["bipartition_sizes"])

    def test_schmidt_entries_are_skipped_above_the_limit(self):
        stats = graph_stats(nx.cycle_graph(8), schmidt_max_nodes=7)
        self.assertIsNone(stats["schmidt_rank_max"])
        self.assertIsNone(stats["schmidt_distribution"])

    def test_orbit_table_row(self):
        g6 = canonical_g6(nx.cycle_graph(6))
        row = orbit_table_row(g6)
        self.assertEqual(len(row), len(ORBIT_TABLE_COLUMNS))
        self.assertEqual(row[0], g6)
        self.assertEqual(row[ORBIT_TABLE_COLUMNS.index("Schmidt_max")], 8)

    def test_orbit_table_rows_in_worker_processes(self):
        graphs = [canonical_g6(nx.cycle_graph(n)) for n in range(3, 9)]
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=2, mp_context=context) as pool:
            rows = list(pool.map(orbit_table_row, graphs))
        self.assertEqual(rows, [orbit_table_row(g6) for g6 in graphs])


class ColouringTests(unittest.TestCase):
    def test_exact_colourability_matches_brute_force(self):
        rng = random.Random(13)
        for _ in range(80):
            graph = nx.gnp_random_graph(rng.randint(1, 7), rng.random(), seed=rng.randrange(10 ** 6))
            for k in (1, 2, 3, 4):
                status, coloring = is_k_colorable(graph, k)
                self.assertEqual(status, brute_force_colorable(graph, k))
                if status:
                    self.assertEqual(set(coloring), set(graph.nodes))
                    self.assertTrue(all(coloring[u] != coloring[v] for u, v in graph.edges()))

    def test_known_graphs_and_step_limit(self):
        self.assertTrue(is_k_colorable(nx.petersen_graph(), 3)[0])
        self.assertFalse(is_k_colorable(nx.petersen_graph(), 2)[0])
        self.assertEqual(is_k_colorable(nx.Graph(), 2), (True, {}))
        self.assertIsNone(is_k_colorable(nx.complete_graph(8), 7, max_steps=3)[0])

    def test_greedy_colouring_is_proper(self):
        for graph in (nx.petersen_graph(), nx.cycle_graph(7), nx.complete_graph(5)):
            coloring = greedy_coloring(graph)
            self.assertTrue(all(coloring[u] != coloring[v] for u, v in graph.edges()))


# ─────────────────────────────────────────────────────────────────────────────
#  GUI (needs a display)
# ─────────────────────────────────────────────────────────────────────────────

def _display_available():
    try:
        import tkinter
        root = tkinter.Tk()
        root.destroy()
        return True
    except Exception:
        return False


@unittest.skipUnless(_display_available(), "no display available")
class GuiTests(unittest.TestCase):
    def setUp(self):
        import tkinter as tk
        from src import graphCreator
        self.tk = tk
        self.gc = graphCreator
        self.errors = []
        self.dialogs = []
        self.root = tk.Tk()
        self.root.geometry("1100x700")
        self.root.report_callback_exception = lambda *exc: self.errors.append(exc)
        for name in ("showinfo", "showwarning", "showerror"):
            patcher = mock.patch(f"tkinter.messagebox.{name}",
                                 lambda *a, _name=name, **k: self.dialogs.append((_name, a)))
            patcher.start()
            self.addCleanup(patcher.stop)
        self.editor = graphCreator.GraphEditor(self.root)
        self.root.update()
        self._last_click = None

    def tearDown(self):
        self.root.destroy()
        self.assertEqual(self.errors, [])

    def run_until(self, condition, timeout=30.0):
        """Run the Tk main loop until ``condition()`` holds or the timeout passes."""
        end = time.time() + timeout

        def tick():
            if condition() or time.time() > end:
                self.root.quit()
            else:
                self.root.after(10, tick)

        self.root.after(10, tick)
        self.root.mainloop()
        return condition()

    def click(self, x, y, widget=None, double=False):
        widget = widget or self.editor.view.canvas
        if not double and self._last_click == (x, y):
            time.sleep(0.7)  # keep separate clicks from forming a double-click
        self._last_click = (x, y)
        for _ in range(2 if double else 1):
            widget.event_generate("<Button-1>", x=x, y=y)
            widget.event_generate("<ButtonRelease-1>", x=x, y=y)
        self.root.update()

    def test_draw_edit_and_local_complementation(self):
        editor = self.editor
        points = [(200, 300), (320, 300), (440, 300), (560, 300)]
        for x, y in points:
            self.click(x, y)
        self.assertEqual([node.label for node in editor.view.nodes], ["0", "1", "2", "3"])
        editor.activate_edge_mode()
        for a, b in [(0, 1), (1, 2), (2, 3)]:
            self.click(*points[a])
            self.click(*points[b])
        self.assertEqual(sorted(editor.view.to_networkx().edges), [(0, 1), (1, 2), (2, 3)])
        self.assertEqual(editor._last_stats["schmidt_distribution"], {2: 5, 4: 2})

        editor.activate_node_mode()
        positions = [(node.x, node.y) for node in editor.view.nodes]
        self.click(*points[1])
        editor.perform_lc()
        expected = graph_from_g6(graphOperations.local_complementation_at_vertex(nx.path_graph(4), 1))
        self.assertEqual(sorted(editor.view.to_networkx().edges), sorted(expected.edges))
        self.assertEqual([(node.x, node.y) for node in editor.view.nodes], positions)

        node = editor.view.nodes[0]
        with mock.patch("tkinter.simpledialog.askstring", return_value="1"):
            self.click(int(node.x), int(node.y), double=True)
        self.assertEqual(node.label, "0")

    def test_load_graph6_and_compute_schmidt(self):
        with mock.patch("tkinter.simpledialog.askstring", return_value=">>graph6<<EFz_\n"):
            self.editor.load_graph6()
        self.assertTrue(nx.is_isomorphic(self.editor.view.to_networkx(),
                                         graph_from_g6("EFz_")))
        self.editor.show_graph(nx.complete_bipartite_graph(3, 3))
        self.editor.compute_schmidt()
        self.assertIn("Maximum Schmidt rank: 4", self.dialogs[-1][1][1])

    def test_orbit_browsing_and_overview(self):
        editor = self.editor
        editor.show_graph(nx.cycle_graph(5))
        editor.compute_orbit()
        self.assertTrue(self.run_until(lambda: not editor._orbit_running))
        expected, _ = reference_orbit(nx.cycle_graph(5))
        self.assertEqual(set(editor.orbit), set(expected))

        visited = []
        for _ in editor.orbit:
            editor.step_orbit(1)
            visited.append(editor.current_g6)
        self.assertEqual(sorted(visited), sorted(editor.orbit))

        editor.show_orbit_overview()
        self.root.update()
        overview = editor.overview
        self.assertEqual(len(overview.view.nodes), len(editor.orbit))
        self.assertEqual([n.label for n in overview.view.nodes if n.marked], [editor.current_g6])
        overview.close()
        editor.show_graph(nx.path_graph(3))
        self.assertIsNone(editor.overview)

    def test_delete_key_in_spinbox_keeps_the_selection(self):
        editor = self.editor
        editor.show_graph(nx.path_graph(3))
        node = editor.view.nodes[0]
        self.click(int(node.x), int(node.y))
        spinbox = next(w for frame in editor.winfo_children() for w in frame.winfo_children()
                       if w.winfo_class() == "TSpinbox")
        spinbox.focus_force()
        self.root.update()
        spinbox.event_generate("<Delete>")
        self.root.update()
        self.assertEqual(editor.view.to_networkx().number_of_nodes(), 3)


if __name__ == "__main__":
    unittest.main()
