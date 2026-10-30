"""Regression tests for the original progressive implementation.

Run from the repository root: python -m unittest discover -s tests -v
Only the existing runtime dependencies are required.
"""
from collections import Counter
import unittest
from unittest.mock import patch

import EdgeAlign as ea


def graph(row, labels, bonds):
    edges = [{"id": f"g{row}e{i}", "type": kind}
             for i, (_, _, kind) in enumerate(bonds)]
    vertices = [{"edges": [edges[i]["id"] for i, (u, v, _) in enumerate(bonds)
                           if vertex in (u, v)], "set_type": label}
                for vertex, label in enumerate(labels)]
    return edges, vertices


def flatten(column):
    if isinstance(column, tuple):
        return [value for part in column for value in flatten(part)]
    return [column]


class CorrectnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ea.configure_logging("CRITICAL")

    def align_and_check(self, graphs, expected_score=None):
        edges, vertices = map(list, zip(*graphs))
        columns, incidence, _, edge_labels, vertex_labels = ea.align_multiple_e_graphs(edges, vertices)
        rows = [flatten(column) for column in columns]
        self.assertTrue(all(len(row) == len(graphs) for row in rows))
        self.assertEqual(len(columns), len(set(columns)))
        self.assertEqual(Counter(incidence), Counter(s for s, _ in vertex_labels))
        self.assertEqual(set(edge_labels), set(map(str, columns)))
        for column in columns:
            self.assertEqual(sum(column in s for s in incidence), 2)
        for i, (source_edges, source_vertices) in enumerate(graphs):
            self.assertEqual(Counter(row[i] for row in rows if row[i] is not None),
                             Counter(e["id"] for e in source_edges))
            projected = Counter()
            for vertex, label in vertex_labels:
                value = frozenset(flatten(c)[i] for c in vertex if flatten(c)[i] is not None)
                if value:
                    projected[(value, label)] += 1
            self.assertEqual(projected, Counter((frozenset(v["edges"]), v["set_type"])
                                                 for v in source_vertices))
        if expected_score is not None:
            self.assertEqual(sum(all(value is not None for value in row) for row in rows),
                             expected_score)

    def test_disconnected_identical_bonds(self):
        self.align_and_check([graph(i, "CCCC", [(0, 1, "SINGLE"), (2, 3, "SINGLE")])
                              for i in range(2)], 2)

    def test_no_compatible_bond_retains_both_inputs(self):
        self.align_and_check([graph(0, "CC", [(0, 1, "SINGLE")]),
                              graph(1, "NN", [(0, 1, "SINGLE")])], 0)

    def test_different_bond_labels_retain_both_inputs(self):
        self.align_and_check([graph(0, "CC", [(0, 1, "SINGLE")]),
                              graph(1, "CC", [(0, 1, "DOUBLE")])], 0)

    def test_isolated_match_vertex_labels_after_refinement(self):
        self.align_and_check([
            graph(0, "CCCNNCC", [(0, 1, "SINGLE"), (1, 2, "SINGLE"),
                  (2, 3, "TRIPLE"), (1, 4, "SINGLE"), (4, 5, "SINGLE"), (4, 6, "SINGLE")]),
            graph(1, "CCCCCCCCO", [(0, 1, "SINGLE"), (1, 2, "AROMATIC"),
                  (2, 3, "AROMATIC"), (3, 4, "AROMATIC"), (4, 5, "AROMATIC"),
                  (5, 6, "AROMATIC"), (6, 7, "SINGLE"), (7, 8, "SINGLE"), (6, 1, "AROMATIC")])], 1)

    def test_unmatched_components_survive_progressive_alignment(self):
        self.align_and_check([graph(i, label * 2, [(0, 1, "SINGLE")])
                              for i, label in enumerate("CNO")], 0)

    def test_implicit_endpoints_both_have_label_entries(self):
        edges = [[{"id": "a", "type": "SINGLE"}], [{"id": "b", "type": "DOUBLE"}]]
        out = ea.align_multiple_e_graphs(edges, [[], []])
        columns, incidence, _, _, labels = out
        self.assertEqual(len(columns), 2)
        self.assertEqual(len(labels), 4)
        for column in columns:
            self.assertEqual(sum(column in vertex for vertex in incidence), 2)

    def test_unexpected_error_is_not_pruned(self):
        inputs = [graph(i, "CC", [(0, 1, "SINGLE")]) for i in range(2)]
        edges, vertices = map(list, zip(*inputs))
        with patch.object(ea, "update_alignment_sets", side_effect=RuntimeError("injected defect")):
            with self.assertRaisesRegex(RuntimeError, "injected defect"):
                ea.align_multiple_e_graphs(edges, vertices)

    def test_expected_conflict_still_prunes_a_candidate(self):
        inputs = [graph(i, "CC", [(0, 1, "SINGLE")]) for i in range(2)]
        edges, vertices = map(list, zip(*inputs))
        with patch.object(ea, "update_alignment_sets", side_effect=ea.AlignmentConflict("test conflict")):
            columns, _, _, _, _ = ea.align_multiple_e_graphs(edges, vertices)
        self.assertEqual(len(columns), 2)
        self.assertTrue(all(None in flatten(column) for column in columns))


if __name__ == "__main__":
    unittest.main()
