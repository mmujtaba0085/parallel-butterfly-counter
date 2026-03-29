"""
Tests for butterfly counting algorithms.

Run with:  pytest test_butterfly.py -v
"""

import pytest
import networkx as nx
from butterfly import butterfly_count, butterfly_count_per_vertex


# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────

def make_bipartite(left, right, edges):
    """Build a bipartite graph with proper node attributes."""
    G = nx.Graph()
    G.add_nodes_from(left, bipartite=0)
    G.add_nodes_from(right, bipartite=1)
    G.add_edges_from(edges)
    return G


# ─────────────────────────────────────────────────────────────
# butterfly_count — correctness
# ─────────────────────────────────────────────────────────────

class TestButterflyCount:

    def test_single_butterfly(self):
        """K_{2,2} contains exactly 1 butterfly."""
        G = make_bipartite([0, 1], [2, 3],
                           [(0, 2), (0, 3), (1, 2), (1, 3)])
        assert butterfly_count(G) == 1

    def test_empty_graph(self):
        G = make_bipartite([0, 1], [2, 3], [])
        assert butterfly_count(G) == 0

    def test_single_edge(self):
        G = make_bipartite([0], [1], [(0, 1)])
        assert butterfly_count(G) == 0

    def test_path_graph(self):
        """0-2-1-3: no butterfly (not K_{2,2})."""
        G = make_bipartite([0, 1], [2, 3],
                           [(0, 2), (1, 2), (1, 3)])
        assert butterfly_count(G) == 0

    def test_k33_has_9_butterflies(self):
        """
        K_{3,3}: C(3,2)^2 / ... let's compute directly.
        Each pair of left nodes share 3 right neighbours → C(3,2)=3 wedges.
        There are C(3,2)=3 left pairs → 3×3 = 9 butterflies.
        """
        G = nx.complete_bipartite_graph(3, 3)
        # nx.complete_bipartite_graph uses 0-based nodes 0..5
        nx.set_node_attributes(G, {i: 0 for i in range(3)}, "bipartite")
        nx.set_node_attributes(G, {i: 1 for i in range(3, 6)}, "bipartite")
        assert butterfly_count(G) == 9

    def test_k44_has_36_butterflies(self):
        """
        K_{4,4}: each pair of left nodes shares 4 right neighbours → C(4,2)=6 butterflies.
        There are C(4,2)=6 left pairs → 6×6 = 36 butterflies.
        """
        G = nx.complete_bipartite_graph(4, 4)
        nx.set_node_attributes(G, {i: 0 for i in range(4)}, "bipartite")
        nx.set_node_attributes(G, {i: 1 for i in range(4, 8)}, "bipartite")
        assert butterfly_count(G) == 36

    def test_two_disjoint_butterflies(self):
        """Two disconnected K_{2,2} components → 2 butterflies."""
        G = nx.Graph()
        G.add_nodes_from([0, 1, 4, 5], bipartite=0)
        G.add_nodes_from([2, 3, 6, 7], bipartite=1)
        G.add_edges_from([(0, 2), (0, 3), (1, 2), (1, 3)])   # butterfly 1
        G.add_edges_from([(4, 6), (4, 7), (5, 6), (5, 7)])   # butterfly 2
        assert butterfly_count(G) == 2

    def test_star_graph_no_butterfly(self):
        """A star K_{1,n} has no butterflies (only one pivot)."""
        G = make_bipartite([0], [1, 2, 3, 4], [(0, 1), (0, 2), (0, 3), (0, 4)])
        assert butterfly_count(G) == 0

    def test_nodes_parameter_both_partitions_give_same_result(self):
        """Pivoting on either partition should give the same count."""
        G = make_bipartite([0, 1], [2, 3],
                           [(0, 2), (0, 3), (1, 2), (1, 3)])
        assert butterfly_count(G, nodes=[0, 1]) == butterfly_count(G, nodes=[2, 3]) == 1

    def test_not_bipartite_raises(self):
        G = nx.cycle_graph(3)   # triangle — not bipartite
        with pytest.raises(nx.NetworkXError):
            butterfly_count(G)

    def test_three_shared_neighbours(self):
        """
        Left: {0,1,2}, Right: {3,4}
        All left nodes connected to both right nodes.
        Each right-node pair (3,4) appears 3 times (once per left node).
        → C(3,2) = 3 butterflies.
        """
        G = make_bipartite([0, 1, 2], [3, 4],
                           [(0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4)])
        assert butterfly_count(G) == 3


# ─────────────────────────────────────────────────────────────
# butterfly_count_per_vertex — correctness
# ─────────────────────────────────────────────────────────────

class TestButterflyCountPerVertex:

    def test_single_butterfly_all_vertices_count_1(self):
        G = make_bipartite([0, 1], [2, 3],
                           [(0, 2), (0, 3), (1, 2), (1, 3)])
        pv = butterfly_count_per_vertex(G)
        assert pv == {0: 1, 1: 1, 2: 1, 3: 1}

    def test_sum_equals_four_times_total(self):
        """Sum of per-vertex counts must equal 4 × total butterfly count."""
        G = nx.complete_bipartite_graph(4, 4)
        nx.set_node_attributes(G, {i: 0 for i in range(4)}, "bipartite")
        nx.set_node_attributes(G, {i: 1 for i in range(4, 8)}, "bipartite")

        total = butterfly_count(G)
        pv = butterfly_count_per_vertex(G)
        assert sum(pv.values()) == 4 * total

    def test_empty_graph_all_zeros(self):
        G = make_bipartite([0, 1], [2, 3], [])
        pv = butterfly_count_per_vertex(G)
        assert all(v == 0 for v in pv.values())

    def test_isolated_node_gets_zero(self):
        G = make_bipartite([0, 1, 99], [2, 3],
                           [(0, 2), (0, 3), (1, 2), (1, 3)])
        pv = butterfly_count_per_vertex(G)
        assert pv[99] == 0

    def test_not_bipartite_raises(self):
        G = nx.cycle_graph(5)
        with pytest.raises(nx.NetworkXError):
            butterfly_count_per_vertex(G)


# ─────────────────────────────────────────────────────────────
# Cross-validate: manual ground truth for K_{3,3}
# ─────────────────────────────────────────────────────────────

class TestK33PerVertex:
    """
    In K_{3,3} every vertex participates in exactly 6 butterflies.

    Proof: pick any vertex u (left side, degree 3).
    For each pair of its 3 right neighbours → 1 butterfly with each other left vertex.
    That's C(3,2)=3 right-pairs × 2 other-left-vertices = 6 butterflies.
    """

    def test_k33_per_vertex(self):
        G = nx.complete_bipartite_graph(3, 3)
        nx.set_node_attributes(G, {i: 0 for i in range(3)}, "bipartite")
        nx.set_node_attributes(G, {i: 1 for i in range(3, 6)}, "bipartite")
        pv = butterfly_count_per_vertex(G)
        assert all(v == 6 for v in pv.values())

    def test_k33_sum(self):
        G = nx.complete_bipartite_graph(3, 3)
        nx.set_node_attributes(G, {i: 0 for i in range(3)}, "bipartite")
        nx.set_node_attributes(G, {i: 1 for i in range(3, 6)}, "bipartite")
        pv = butterfly_count_per_vertex(G)
        # 9 butterflies × 4 vertices each = 36
        assert sum(pv.values()) == 36
