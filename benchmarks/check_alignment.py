"""Independent checks used in the original Appendix B measurements.

Extracted unchanged from the archived study utility; no RDKit dependency.
These checks run after the timed alignment call.
"""
from collections import Counter


def flatten(column):
    if isinstance(column, (tuple, list)):
        return [v for x in column for v in flatten(x)]
    return [column]

def validate_output(output, inputs):
    columns, sets, _, edge_types, set_types = output
    source_edges, source_sets = inputs
    n = len(source_edges)
    rows = [flatten(c) for c in columns]
    errors = []
    if any(len(c) != n for c in rows):
        return {"valid": False, "errors": ["wrong number of rows in a column"]}
    if len(set(columns)) != len(columns):
        errors.append("duplicate alignment columns")
    if any(all(v is None for v in c) for c in rows):
        errors.append("gap-only column")
    for i in range(n):
        observed = Counter(c[i] for c in rows if c[i] is not None)
        expected = Counter(e["id"] for e in source_edges[i])
        if observed != expected:
            errors.append(f"row {i}: bonds not recovered exactly")
    for c in columns:
        if sum(c in s for s in sets) != 2:
            errors.append("column does not have exactly two endpoints")
            break
    for i, a in enumerate(sets):
        if not a or not set(a).issubset(columns):
            errors.append("empty vertex or unknown column")
        if any(len(a & b) > 1 for b in sets[i + 1:]):
            errors.append("parallel edges in the alignment graph")
            break
    if Counter(sets) != Counter(s for s, _ in set_types):
        errors.append("vertex labels do not cover the alignment vertices")
    if set(edge_types) != set(map(str, columns)):
        errors.append("bond labels do not cover the alignment columns")
    for i in range(n):
        expected = Counter((frozenset(s["edges"]), s["set_type"]) for s in source_sets[i])
        observed = Counter()
        edge_labels = {e["id"]: e["type"] for e in source_edges[i]}
        for s, typ in set_types:
            projected = frozenset(flatten(c)[i] for c in s if flatten(c)[i] is not None)
            if projected:
                observed[(projected, typ)] += 1
        if observed != expected:
            errors.append(f"row {i}: labeled vertices not recovered exactly")
        for c, row in zip(columns, rows):
            if row[i] is not None and edge_types.get(str(c)) != edge_labels[row[i]]:
                errors.append(f"row {i}: bond label mismatch")
                break
    occupancy = [sum(x is not None for x in c) for c in rows]
    return {"valid": not errors, "errors": sorted(set(errors)), "columns": len(columns),
            "vertices": len(sets), "core_bonds": sum(k == n for k in occupancy),
            "sum_pairs": sum(k * (k - 1) // 2 for k in occupancy),
            "occupancy": occupancy}
