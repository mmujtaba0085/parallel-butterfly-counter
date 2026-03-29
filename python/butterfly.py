"""
Butterfly counting algorithms for bipartite graphs.

A butterfly is a complete bipartite subgraph K_{2,2} — four vertices
(two from each partition) all connected to each other, forming a
"cycle of length 4" or equivalently a pair of wedges sharing both endpoints.

This module provides:
- butterfly_count(G)             : total butterfly count in the graph
- butterfly_count_per_vertex(G)  : per-vertex butterfly participation count

Algorithm
---------
Based on wedge-based intersection (Sanei-Mehri et al., KDD 2018):

1. Rank vertices by degree (ascending) — processes low-degree vertices first
   to minimise redundant wedge enumeration.
2. For each vertex v in the ranked order, enumerate all pairs of neighbours
   (u1, u2) — these form a "wedge" centred at v.
3. Count how many times each (u1, u2) pair appears. If a pair appears k times,
   it contributes C(k, 2) = k*(k-1)//2 butterflies.

Time complexity : O(sum_{v} d(v)^2) where d(v) is the degree of v.
Space complexity: O(E) for the wedge map.

References
----------
.. [1] Sanei-Mehri, S. V., Sariyuce, A. E., & Tirthapura, S. (2018).
       Butterfly counting in bipartite networks.
       In Proceedings of the 24th ACM SIGKDD (pp. 2150-2159).
       https://doi.org/10.1145/3219819.3220097
"""

from collections import defaultdict
from itertools import combinations

import networkx as nx
from networkx.algorithms.bipartite import sets as bipartite_sets

__all__ = ["butterfly_count", "butterfly_count_per_vertex"]


def _check_bipartite(G):
    """Raise if G is not bipartite."""
    if not nx.is_bipartite(G):
        raise nx.NetworkXError(
            "butterfly_count requires a bipartite graph. "
            "Use nx.is_bipartite(G) to check."
        )


def _get_partition(G, nodes):
    """Return a node set for the pivot partition.

    Priority order:
    1. *nodes* argument (user-supplied)
    2. ``bipartite`` node attribute (0 → left partition)
    3. 2-coloring via BFS (handles disconnected graphs correctly)
    """
    if nodes is not None:
        return set(nodes)

    # Try node attribute first (standard networkx bipartite convention)
    attr = nx.get_node_attributes(G, "bipartite")
    if attr:
        left = {n for n, v in attr.items() if v == 0}
        right = {n for n, v in attr.items() if v == 1}
        # Pivot on the smaller side
        return left if len(left) <= len(right) else right

    # Fall back to BFS 2-coloring — works on disconnected graphs
    color = {}
    for start in G.nodes():
        if start in color:
            continue
        queue = [start]
        color[start] = 0
        while queue:
            v = queue.pop()
            for nbr in G.neighbors(v):
                if nbr not in color:
                    color[nbr] = 1 - color[v]
                    queue.append(nbr)

    left = {n for n, c in color.items() if c == 0}
    right = {n for n, c in color.items() if c == 1}
    return left if len(left) <= len(right) else right


def _rank_vertices(G, nodes):
    """Return *nodes* sorted by ascending degree (ties broken by node id).

    Processing low-degree vertices first keeps the wedge map small and
    matches the vertex-ranking strategy in the PARBUTTERFLY paper.
    """
    return sorted(nodes, key=lambda v: (G.degree(v), v))


def butterfly_count(G, nodes=None):
    """Count the total number of butterflies (K_{2,2} subgraphs) in *G*.

    Parameters
    ----------
    G : NetworkX Graph
        An undirected bipartite graph.  The graph must have a ``bipartite``
        node attribute (0 or 1) as produced by the bipartite graph generators
        or set manually.
    nodes : container, optional
        Nodes from one partition of *G* to use as the "pivot" side for wedge
        enumeration.  If *None*, the partition with fewer nodes is chosen
        automatically (usually faster).

    Returns
    -------
    int
        Total number of butterflies.  Each butterfly is counted once.

    Raises
    ------
    NetworkXError
        If *G* is not bipartite.

    Examples
    --------
    Build a simple butterfly (K_{2,2}):

    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_nodes_from([0, 1], bipartite=0)   # left partition
    >>> G.add_nodes_from([2, 3], bipartite=1)   # right partition
    >>> G.add_edges_from([(0, 2), (0, 3), (1, 2), (1, 3)])
    >>> from butterfly import butterfly_count
    >>> butterfly_count(G)
    1

    A graph with two overlapping butterflies:

    >>> G2 = nx.complete_bipartite_graph(3, 3)
    >>> butterfly_count(G2)
    9

    Notes
    -----
    The algorithm centres wedge enumeration on one partition (the pivot side).
    For each pivot vertex *v*, every pair of neighbours (u1, u2) in the
    *opposite* partition forms a wedge.  When two different pivots share the
    same pair (u1, u2) the four vertices {v, v', u1, u2} form a butterfly.

    The formula C(k, 2) converts a wedge-pair count *k* to a butterfly count
    in O(1) per pair, giving overall O(sum d(v)^2) complexity.
    """
    _check_bipartite(G)

    if G.number_of_edges() == 0:
        return 0

    pivot_nodes = _get_partition(G, nodes)

    # Rank pivot vertices by ascending degree
    ranked = _rank_vertices(G, pivot_nodes)

    # wedge_counts[(u1, u2)] = number of pivot vertices connected to both u1 & u2
    wedge_counts = defaultdict(int)

    for v in ranked:
        # Neighbours of v are in the OPPOSITE partition
        nbrs = list(G.neighbors(v))
        for u1, u2 in combinations(nbrs, 2):
            # Canonical ordering so (u1,u2) and (u2,u1) map to same key
            key = (u1, u2) if u1 < u2 else (u2, u1)
            wedge_counts[key] += 1

    # Each wedge-pair count k → C(k,2) butterflies
    total = 0
    for k in wedge_counts.values():
        if k >= 2:
            total += k * (k - 1) // 2

    return total


def butterfly_count_per_vertex(G, nodes=None):
    """Count how many butterflies each vertex *participates in*.

    A vertex participates in a butterfly for every K_{2,2} subgraph that
    contains it.  The sum of all per-vertex counts equals 4× the total
    butterfly count (each butterfly has four vertices).

    Parameters
    ----------
    G : NetworkX Graph
        An undirected bipartite graph.
    nodes : container, optional
        Pivot partition for wedge enumeration (same semantics as
        :func:`butterfly_count`).

    Returns
    -------
    dict
        Mapping ``{node: count}`` for every node in *G*.

    Raises
    ------
    NetworkXError
        If *G* is not bipartite.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_nodes_from([0, 1], bipartite=0)
    >>> G.add_nodes_from([2, 3], bipartite=1)
    >>> G.add_edges_from([(0, 2), (0, 3), (1, 2), (1, 3)])
    >>> from butterfly import butterfly_count_per_vertex
    >>> butterfly_count_per_vertex(G)
    {0: 1, 1: 1, 2: 1, 3: 1}
    """
    _check_bipartite(G)

    per_vertex = defaultdict(int)

    if G.number_of_edges() == 0:
        return dict(per_vertex)

    pivot_nodes = _get_partition(G, nodes)
    ranked = _rank_vertices(G, pivot_nodes)

    wedge_counts = defaultdict(int)
    # Also track which pivots contributed to each wedge pair
    wedge_pivots = defaultdict(list)

    for v in ranked:
        nbrs = list(G.neighbors(v))
        for u1, u2 in combinations(nbrs, 2):
            key = (u1, u2) if u1 < u2 else (u2, u1)
            wedge_counts[key] += 1
            wedge_pivots[key].append(v)

    for (u1, u2), k in wedge_counts.items():
        if k >= 2:
            # Number of butterflies this pair (u1,u2) participates in
            bf = k * (k - 1) // 2
            # u1 and u2 are each in bf butterflies
            per_vertex[u1] += bf
            per_vertex[u2] += bf
            # Each pivot vertex that contributed a wedge is in (k-1) butterflies
            # (it pairs with every OTHER pivot that shares this endpoint pair)
            for v in wedge_pivots[(u1, u2)]:
                per_vertex[v] += k - 1

    # Ensure every node has an entry (even isolated ones)
    for node in G.nodes():
        if node not in per_vertex:
            per_vertex[node] = 0

    return dict(per_vertex)
