"""Local complementation and LC-orbit exploration for graph states.

Local complementation (LC) at a vertex v replaces the subgraph induced on the
neighbourhood N(v) by its complement. On graph states, LC at v is implemented
by a local Clifford unitary, and two graph states are local-Clifford
equivalent if and only if their graphs are related by a sequence of LC
operations [Van den Nest, Dehaene & De Moor, Phys. Rev. A 69, 022316 (2004)].

Orbits are explored on canonical graph6 strings (see ``toolset``), so graphs
are identified up to isomorphism. An orbit returned here is therefore the set
of isomorphism classes reachable by LC, i.e. the graph states equivalent to
the seed under local Clifford operations combined with qubit relabelling.
"""

from collections import deque
import operator
import threading
import warnings

import networkx as nx
import numpy as np
from tqdm import TqdmExperimentalWarning
from tqdm.tk import tqdm_tk

from src.toolset import canonical_relabel_n as canonical_relabel


_GRAPH6_HEADER = ">>graph6<<"


# ─────────────────────────────────────────────────────────────────────────────
#  Progress reporting
# ─────────────────────────────────────────────────────────────────────────────

class _TkProgressProxy:
    """Run every call on a ``tqdm_tk`` progress bar in the Tk main thread.

    Tk widgets may only be touched from the thread that runs the mainloop.
    Orbit exploration runs in a worker thread, so each call is scheduled with
    ``master.after`` and the worker waits until it has been executed. Without
    a ``master``, or when already on the main thread, calls run directly.
    """

    def __init__(self, master, **kwargs):
        self._master = master
        self._pbar = None
        self._call(lambda: self._create(kwargs))

    def _create(self, kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", TqdmExperimentalWarning)
            self._pbar = tqdm_tk(**kwargs)

    def _should_dispatch(self):
        return (
            self._master is not None
            and threading.current_thread() is not threading.main_thread()
        )

    def _call(self, callback):
        if not self._should_dispatch():
            return callback()

        done = threading.Event()
        result = []
        error = []

        def run():
            try:
                result.append(callback())
            except BaseException as exc:
                error.append(exc)
            finally:
                done.set()

        self._master.after(0, run)
        done.wait()
        if error:
            raise error[0]
        return result[0] if result else None

    @property
    def n(self):
        return self._pbar.n

    @property
    def total(self):
        return self._pbar.total

    @total.setter
    def total(self, total):
        self._call(lambda: self._set_total(total))

    def _set_total(self, total):
        self._pbar.total = total
        tk_pbar = getattr(self._pbar, "_tk_pbar", None)
        if tk_pbar is not None:
            tk_pbar.configure(maximum=total)

    def set_description(self, text):
        self._call(lambda: self._pbar.set_description(text))

    def update(self, n=1):
        self._call(lambda: self._pbar.update(n))

    def finish(self, description):
        """Mark the bar as complete and show a final description."""
        self._call(lambda: self._finish(description))

    def _finish(self, description):
        if self._pbar.n == 0:
            self._set_total(1)
            self._pbar.update(1)
        elif self._pbar.total != self._pbar.n:
            self._set_total(self._pbar.n)
        self._pbar.set_description(description)
        self._pbar.refresh()

    def close(self):
        self._call(self._pbar.close)


class _SilentProgress:
    """Progress counter with the interface of _TkProgressProxy and no window."""

    def __init__(self, total):
        self.n = 0
        self.total = total

    def set_description(self, text):
        pass

    def update(self, n=1):
        self.n += n

    def finish(self, description):
        pass

    def close(self):
        pass


# ─────────────────────────────────────────────────────────────────────────────
#  Local complementation
# ─────────────────────────────────────────────────────────────────────────────

def _adjacency(graph):
    """Return the 0/1 adjacency matrix of ``graph`` in ``graph.nodes`` order.

    Edge attributes such as ``weight`` are ignored.
    """
    if graph.is_directed() or graph.is_multigraph():
        raise TypeError("expected a simple undirected nx.Graph")
    if nx.number_of_selfloops(graph) > 0:
        raise ValueError("graph states are defined on graphs without self-loops")
    return nx.to_numpy_array(graph, dtype=np.uint8, weight=None)


def _complement_neighbourhood(A, v):
    """Return a copy of ``A`` with LC applied at index ``v``.

    Returns ``None`` when ``v`` has fewer than two neighbours, because LC
    then leaves the graph unchanged.
    """
    nbrs = np.flatnonzero(A[v])
    if nbrs.size < 2:
        return None
    M = A.copy()
    block = np.ix_(nbrs, nbrs)
    M[block] ^= 1
    np.fill_diagonal(M, 0)
    return M


def _matrix_to_graph6(M):
    """graph6 string of the graph with adjacency matrix ``M``, vertex order kept."""
    g6 = nx.to_graph6_bytes(nx.from_numpy_array(M), header=False)
    return g6.decode("ascii").rstrip("\n")


def _parse_graph6(g6):
    """Read a graph6 string or bytes object, with or without header/newline."""
    if isinstance(g6, (bytes, bytearray)):
        g6 = bytes(g6).decode("ascii")
    g6 = g6.strip()
    if g6.startswith(_GRAPH6_HEADER):
        g6 = g6[len(_GRAPH6_HEADER):]
    return nx.from_graph6_bytes(g6.encode("ascii"))


def local_complementation_at_vertex(graph, v):
    """Apply LC at a single vertex and return the result as graph6.

    Parameters
    ----------
    graph : networkx.Graph
        Simple undirected graph.
    v : int
        Position of the vertex in ``list(graph.nodes)`` (0-based), not its label.

    Returns
    -------
    str
        graph6 string (no header, no newline). Vertex ``i`` of the result is
        the ``i``-th vertex of ``graph.nodes``. The result is not canonicalised,
        so a drawing of the graph can keep its vertex positions.

    Raises
    ------
    IndexError
        If ``v`` is not a valid vertex position.
    """
    A = _adjacency(graph)
    n = A.shape[0]
    v = operator.index(v)
    if not 0 <= v < n:
        raise IndexError(f"vertex position {v} is out of range for {n} vertices")
    M = _complement_neighbourhood(A, v)
    return _matrix_to_graph6(A if M is None else M)


def localComplementationOperationUI(graph):
    """Apply LC at every vertex of ``graph`` and collect the distinct results.

    Vertices are visited in ``graph.nodes`` order. Vertices with fewer than
    two neighbours are skipped. Results are canonicalised, so isomorphic
    outcomes are reported once. A result can be isomorphic to ``graph`` itself.

    Returns
    -------
    targets : list of str
        Canonical graph6 strings of the distinct LC neighbours, in the order
        they were found.
    operations : list of (str, int)
        Pairs ``(target, v)``, where ``v`` is the position in ``graph.nodes``
        of the first vertex whose LC produces ``target``.
    """
    A = _adjacency(graph)
    seen = set()
    targets = []
    operations = []
    for v in range(A.shape[0]):
        M = _complement_neighbourhood(A, v)
        if M is None:
            continue
        target = canonical_relabel(nx.from_numpy_array(M)).decode("ascii")
        if target not in seen:
            seen.add(target)
            targets.append(target)
            operations.append((target, v))
    return targets, operations


def localComplementationOperation(graph):
    """Canonical graph6 strings of the distinct LC neighbours of ``graph``.

    Same as :func:`localComplementationOperationUI` without the vertex record.
    """
    targets, _operations = localComplementationOperationUI(graph)
    return targets


# ─────────────────────────────────────────────────────────────────────────────
#  Orbit exploration
# ─────────────────────────────────────────────────────────────────────────────

def exploreOrbit(seeds, master=None, progress=True):
    """Explore the LC orbit of one or more seed graphs by breadth-first search.

    Every graph is stored by its canonical graph6 string, so the search runs
    over isomorphism classes and terminates once no new class appears.

    Parameters
    ----------
    seeds : iterable of str or bytes
        graph6 encodings of the seed graphs. A ``>>graph6<<`` header and a
        trailing newline are accepted. The input is not modified.
    master : tkinter widget or None
        Tk root that owns the progress window. Pass the root when calling from
        a worker thread of the GUI.
    progress : bool
        Show a Tk progress window. Use ``False`` in scripts and on machines
        without a display.

    Returns
    -------
    orbit : list of str
        Canonical graph6 strings of all graphs reached, in BFS order. The
        first entry is the canonical form of the first seed. If the seeds lie
        in different orbits, the result is the union of their orbits.
    neighbors : dict
        Maps every string in ``orbit`` to its list of LC moves ``(target, v)``
        as returned by :func:`localComplementationOperationUI`. ``v`` is a
        vertex of the canonical graph given by the key.
    """
    start = [canonical_relabel(_parse_graph6(s)).decode("ascii") for s in seeds]

    if progress:
        pbar = _TkProgressProxy(
            master,
            total=max(len(start), 1),
            desc="Exploring LC orbit",
            unit="graph",
            tk_parent=master,
            leave=False,  # close the window when the search is finished
        )
    else:
        pbar = _SilentProgress(total=max(len(start), 1))
    orbit = []
    neighbors = {}
    completed = False
    try:
        for seed in start:
            if seed in neighbors:
                if pbar.n < (pbar.total or 0):
                    pbar.update(1)
                continue

            queue = deque([seed])
            queued = {seed}
            while queue:
                current = queue.popleft()
                queued.discard(current)

                targets, operations = localComplementationOperationUI(
                    _parse_graph6(current)
                )
                orbit.append(current)
                neighbors[current] = operations

                for target in targets:
                    if target not in neighbors and target not in queued:
                        queue.append(target)
                        queued.add(target)

                pending = pbar.n + 1 + len(queue)
                if pending > (pbar.total or 0):
                    pbar.total = pending
                pbar.set_description(
                    f"Orbit size: {len(orbit)} | Queue: {len(queue)}"
                )
                pbar.update(1)
        completed = True
    finally:
        if completed:
            pbar.finish(f"Orbit complete: {len(orbit)}")
        pbar.close()

    return orbit, neighbors
