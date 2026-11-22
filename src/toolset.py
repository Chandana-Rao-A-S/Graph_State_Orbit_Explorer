"""Canonical labelling of simple undirected graphs.

The LC-orbit explorer stores every graph as the graph6 string of its
canonical form. Two graphs are then the same orbit member exactly when
their strings are equal. Canonical forms are computed by nauty through
the pynauty bindings.
"""

import networkx as nx
import pynauty as pyn


def canonical_relabel_n(nx_g: nx.Graph) -> bytes:
    """Return the graph6 encoding of the canonical form of ``nx_g``.

    Two graphs receive the same string if and only if they are isomorphic.
    Node labels may be of any hashable type (integers, strings, gaps in the
    numbering); only the graph structure is used. Edge attributes are ignored.

    Parameters
    ----------
    nx_g : networkx.Graph
        Simple undirected graph without self-loops.

    Returns
    -------
    bytes
        graph6 string without the ``>>graph6<<`` header and without a
        trailing newline.

    Raises
    ------
    TypeError
        If ``nx_g`` is directed or a multigraph.
    ValueError
        If ``nx_g`` contains a self-loop.
    """
    if nx_g.is_directed() or nx_g.is_multigraph():
        raise TypeError("canonical_relabel_n expects a simple undirected nx.Graph")
    if nx.number_of_selfloops(nx_g) > 0:
        raise ValueError("canonical_relabel_n expects a graph without self-loops")

    n = nx_g.number_of_nodes()
    if n == 0:
        # pynauty cannot build a graph with zero vertices.
        return nx.to_graph6_bytes(nx.empty_graph(0), header=False).rstrip(b"\n")

    # The canonical form depends only on the isomorphism class, so mapping
    # arbitrary labels onto 0..n-1 first cannot change the result.
    g = nx.convert_node_labels_to_integers(nx_g)
    adjacency = {v: list(g.neighbors(v)) for v in g.nodes}
    canon = pyn.canon_graph(
        pyn.Graph(number_of_vertices=n, directed=False, adjacency_dict=adjacency)
    )

    # Building an undirected nx.Graph makes the result symmetric regardless of
    # whether pynauty lists each edge once or twice.
    canonical = nx.empty_graph(n)
    for u, nbrs in canon.adjacency_dict.items():
        canonical.add_edges_from((u, v) for v in nbrs)

    return nx.to_graph6_bytes(canonical, header=False).rstrip(b"\n")
