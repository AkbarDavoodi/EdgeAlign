import argparse
import logging
import json
from collections import defaultdict, deque, Counter
from itertools import combinations
from copy import deepcopy
import networkx as nx
import matplotlib.pyplot as plt
import itertools
from typing import Any, List, Set, Tuple, Dict, Optional, FrozenSet


def configure_logging(log_level: str, log_file: Optional[str] = None) -> None:
    """Configures logging based on command-line arguments."""
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=numeric_level,
        format="%(levelname)s: %(message)s",
        handlers=handlers
    )
    # Create a logger for our module
    global logger
    logger = logging.getLogger(__name__)


def pad_none(depth: int) -> Optional[Any]:
    if depth <= 1:
        return None
    else:
        return (pad_none(depth - 1), None)


def pad_edge(e: Any, current_depth: int) -> Any:
    if current_depth <= 1:
        return e
    else:
        return (pad_none(current_depth - 1), e)


def refine_alignment_with_unaligned_edges(
        alignment_edges: List[Tuple[Any, Any]],
        alignment_sets: List[FrozenSet[Tuple[Any, Any]]],
        set_matches: Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]],
        sets_1: List[Set[Any]],
        sets_2: List[Set[Any]],
        edges_1: Set[Any],
        edges_2: Set[Any],
        current_depth: int
) -> Tuple[List[Tuple[Any, Any]], List[FrozenSet[Tuple[Any, Any]]], Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]]]:
    logger.debug("=== START PROCESSING UNALIGNED EDGES IN FIRST E-GRAPH ===")
    logger.debug(f"set_matches: {set_matches}")
    logger.debug(f"alignment_edges: {alignment_edges}")

    aligned_edges_1 = {e1 for e1, _ in alignment_edges}
    U1 = set(edges_1) - aligned_edges_1

    logger.debug(f"Initial unaligned edges in first E-graph (U1): {U1}")

    while U1:
        logger.debug(f"U1 = {U1}")
        progress = False

        for e in list(U1):
            endpoint_sets = [s for s in sets_1 if e in s]
            logger.debug(f"endpoint_sets = {endpoint_sets}")
            logger.debug(f"Current set_matches: {set_matches}")

            if len(endpoint_sets) != 2:
                logger.error(f"Edge {e} should appear in exactly two sets, but found {len(endpoint_sets)}")
                continue

            endpoint_1, endpoint_2 = endpoint_sets

            w1 = next((key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_1)), None)
            w2 = next((key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_2)), None)

            logger.debug(f"w1 = {w1}\nw2 = {w2}")

            if w1 and w2:
                logger.debug(f"Adding ({e}, None) to two related sets.")
                updated_w1, updated_w2 = set(w1), set(w2)
                updated_w1.add((e, None))
                updated_w2.add((e, None))
                new_w1, new_w2 = frozenset(updated_w1), frozenset(updated_w2)
                alignment_sets[alignment_sets.index(w1)] = new_w1
                alignment_sets[alignment_sets.index(w2)] = new_w2

                matching_entry_1 = next((key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_1)),
                                        None)
                matching_entry_2 = next((key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_2)),
                                        None)

                if matching_entry_1:
                    old_value = set_matches[matching_entry_1]
                    del set_matches[matching_entry_1]
                    set_matches[new_w1] = old_value
                    logger.debug(f"Updated set_matches: Removed {matching_entry_1}, Added {new_w1} -> {old_value}")
                else:
                    logger.error(f"No matching set found for endpoint_1: {endpoint_1}")

                if matching_entry_2:
                    old_value = set_matches[matching_entry_2]
                    del set_matches[matching_entry_2]
                    set_matches[new_w2] = old_value
                    logger.debug(f"Updated set_matches: Removed {matching_entry_2}, Added {new_w2} -> {old_value}")
                else:
                    logger.error(f"No matching set found for endpoint_2: {endpoint_2}")

                alignment_edges.append((e, None))
                U1.remove(e)
                progress = True

            elif w1 or w2:
                logger.debug(f"Adding ({e}, None) to one related set and creating a new set.")
                existing_w = w1 if w1 else w2
                updated_w = set(existing_w)
                updated_w.add((e, None))
                new_w = frozenset(updated_w)
                alignment_sets[alignment_sets.index(existing_w)] = new_w

                matching_entry = next(
                    (key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_1 if w1 else endpoint_2)),
                    None
                )
                logger.debug(f"matching_entry = {matching_entry}")
                if matching_entry:
                    old_value = set_matches[matching_entry]
                    del set_matches[matching_entry]
                    set_matches[new_w] = old_value
                    logger.debug(f"Updated set_matches: Removed {matching_entry}, Added {new_w} -> {old_value}")
                else:
                    logger.error("No matching set found for endpoint_1 or endpoint_2.")

                new_other_w = frozenset({(e, None)})
                alignment_sets.append(new_other_w)

                if w1:
                    set_matches[new_other_w] = (endpoint_2, "nothing")
                    logger.debug(f"Created new set_matches entry: {new_other_w} -> ({endpoint_2}, 'nothing')")
                else:
                    set_matches[new_other_w] = (endpoint_1, "nothing")
                    logger.debug(f"Created new set_matches entry: {new_other_w} -> ({endpoint_1}, 'nothing')")

                alignment_edges.append((e, None))
                U1.remove(e)
                progress = True

            else:
                logger.debug(f"Could not find a related set for {e} in first E-graph. It is postponed.")
                continue

        if not progress:
            logger.error(f"Could not find a related set for any of {U1} in first E-graph.")
            break

    logger.debug("=== START PROCESSING UNALIGNED EDGES IN SECOND E-GRAPH ===")
    logger.debug(f"alignment_sets = {alignment_sets}")
    logger.debug(f"alignment_edges = {alignment_edges}")
    logger.debug(f"set_matches = {set_matches}")

    aligned_edges_2 = {e2 for _, e2 in alignment_edges}
    U2 = set(edges_2) - aligned_edges_2

    logger.debug(f"Initial unaligned edges in second E-graph (U2): {U2}")

    while U2:
        logger.debug(f"U2 = {U2}")
        progress = False

        for e in list(U2):
            endpoint_sets = [s for s in sets_2 if e in s]
            logger.debug(f"endpoint_sets = {endpoint_sets}")
            logger.debug(f"Current set_matches: {set_matches}")

            if len(endpoint_sets) != 2:
                logger.error(f"Edge {e} should appear in exactly two sets, but found {len(endpoint_sets)}")
                continue

            endpoint_1, endpoint_2 = endpoint_sets

            w1 = next((key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_1)), None)
            w2 = next((key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_2)), None)
            logger.debug(f"w1 = {w1}\nw2 = {w2}")

            if w1 and w2:
                logger.debug(f"Adding (None, {e}) to two related sets.")
                updated_w1, updated_w2 = set(w1), set(w2)
                updated_w1.add(pad_edge(e, current_depth))
                updated_w2.add(pad_edge(e, current_depth))
                new_w1, new_w2 = frozenset(updated_w1), frozenset(updated_w2)
                alignment_sets[alignment_sets.index(w1)] = new_w1
                alignment_sets[alignment_sets.index(w2)] = new_w2

                matching_entry_1 = next((key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_1)),
                                        None)
                matching_entry_2 = next((key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_2)),
                                        None)
                if matching_entry_1:
                    old_value = set_matches[matching_entry_1]
                    del set_matches[matching_entry_1]
                    set_matches[new_w1] = old_value
                    logger.debug(f"Updated set_matches: Removed {matching_entry_1}, Added {new_w1} -> {old_value}")
                else:
                    logger.error(f"No matching set found for endpoint_1: {endpoint_1}")

                if matching_entry_2:
                    old_value = set_matches[matching_entry_2]
                    del set_matches[matching_entry_2]
                    set_matches[new_w2] = old_value
                    logger.debug(f"Updated set_matches: Removed {matching_entry_2}, Added {new_w2} -> {old_value}")
                else:
                    logger.error(f"No matching set found for endpoint_2: {endpoint_2}")

                alignment_edges.append(pad_edge(e, current_depth))
                U2.remove(e)
                progress = True

            elif w1 or w2:
                logger.debug(f"Adding (None, {e}) to one related set and creating a new set.")
                existing_w = w1 if w1 else w2
                updated_w = set(existing_w)
                updated_w.add(pad_edge(e, current_depth))
                new_w = frozenset(updated_w)
                alignment_sets[alignment_sets.index(existing_w)] = new_w

                matching_entry = next(
                    (key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_1 if w1 else endpoint_2)),
                    None
                )
                logger.debug(f"matching_entry = {matching_entry}")
                if matching_entry:
                    old_value = set_matches[matching_entry]
                    del set_matches[matching_entry]
                    set_matches[new_w] = old_value
                    logger.debug(f"Updated set_matches: Removed {matching_entry}, Added {new_w} -> {old_value}")
                else:
                    logger.error("No matching set found for endpoint_2 or endpoint_1.")

                new_other_w = frozenset({pad_edge(e, current_depth)})
                alignment_sets.append(new_other_w)

                if w1:
                    set_matches[new_other_w] = ("nothing", endpoint_2)
                    logger.debug(f"Created new set_matches entry: {new_other_w} -> ('nothing', {endpoint_2})")
                else:
                    set_matches[new_other_w] = ("nothing", endpoint_1)
                    logger.debug(f"Created new set_matches entry: {new_other_w} -> ('nothing', {endpoint_1})")

                alignment_edges.append(pad_edge(e, current_depth))
                U2.remove(e)
                progress = True

            else:
                logger.debug(f"Could not find a related set for {e} in second E-graph. It is postponed.")
                continue

        if not progress:
            logger.error(f"Could not find a related set for any of {U2} in second E-graph.")
            break

    logger.debug("=== FINISHED PROCESSING UNALIGNED EDGES ===")
    logger.debug(f"Final alignment_edges: {alignment_edges}")
    logger.debug(f"Final alignment_sets: {alignment_sets}")
    logger.debug(f"Final set_matches: {set_matches}")

    return alignment_edges, alignment_sets, set_matches


def check_implied_edge_matches(
        set_matches: Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]],
        alignment_edges: List[Tuple[Any, Any]],
        alignment_sets: List[FrozenSet[Tuple[Any, Any]]]
) -> Tuple[Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]], List[Tuple[Any, Any]], List[FrozenSet[Tuple[Any, Any]]]]:
    logger.debug("Checking for implied edge matches...")
    updated = True
    while updated:
        updated = False
        for (w1, (v1, u1)), (w2, (v2, u2)) in combinations(set_matches.items(), 2):
            if w1 != w2:
                v1_set, u1_set = frozenset(v1), frozenset(u1)
                v2_set, u2_set = frozenset(v2), frozenset(u2)
                if v1_set.intersection(v2_set) and u1_set.intersection(u2_set):
                    e1 = list(v1_set.intersection(v2_set))[0]
                    e2 = list(u1_set.intersection(u2_set))[0]
                    if not any((e1, e2) in s for s in alignment_sets):
                        logger.debug(f"Found implied edge: ({e1}, {e2})")
                        logger.debug(f"set_matches = {set_matches}")
                        alignment_edges.append((e1, e2))
                        updated_w1 = set(w1)
                        updated_w1.add((e1, e2))
                        new_w1 = frozenset(updated_w1)
                        updated_w2 = set(w2)
                        updated_w2.add((e1, e2))
                        new_w2 = frozenset(updated_w2)
                        alignment_sets[alignment_sets.index(w1)] = new_w1
                        alignment_sets[alignment_sets.index(w2)] = new_w2
                        if frozenset(w1) in set_matches:
                            set_matches[frozenset(new_w1)] = set_matches.pop(frozenset(w1))
                        if frozenset(w2) in set_matches:
                            set_matches[frozenset(new_w2)] = set_matches.pop(frozenset(w2))
                        logger.debug(f"Updated w1: {new_w1}")
                        logger.debug(f"Updated w2: {new_w2}")
                        updated = True
                        break
        if updated:
            logger.debug("One round complete. Restarting iteration due to updates.")
        else:
            logger.debug("No new implied edge found in this round. Terminating.")
    logger.debug(f"Final alignment_edges: {alignment_edges}")
    logger.debug(f"Final alignment_sets: {alignment_sets}")
    logger.debug(f"Final set_matches: {set_matches}")
    return set_matches, alignment_edges, alignment_sets


def update_alignment_sets(
        e1: Any,
        e2: Any,
        alignment_sets: List[FrozenSet[Tuple[Any, Any]]],
        set_matches: Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]],
        sets_1: List[Set[Any]],
        sets_2: List[Set[Any]],
        isolated_edges: Dict[Tuple[Any, Any], Tuple[Any, Any, Any, Any]],
        current_alignment: List[Tuple[Any, Any]]
) -> Tuple[List[FrozenSet[Tuple[Any, Any]]], Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]]]:
    logger.debug(f"Updating alignment sets for edge ({e1}, {e2})...")
    f1, f2 = None, None
    is_isolated = True
    for alignment_set in alignment_sets:
        for pair in alignment_set:
            if any(e1 in s and pair[0] in s for s in sets_1) and any(e2 in s and pair[1] in s for s in sets_2):
                f1, f2 = pair
                is_isolated = False
                break
        if not is_isolated:
            break
    if is_isolated:
        logger.debug("Edge is isolated. Deferring decision.")
        v1_candidates = [s for s in sets_1 if e1 in s]
        u1_candidates = [s for s in sets_2 if e2 in s]
        v1 = frozenset(v1_candidates[0]) if v1_candidates else frozenset()
        v2 = frozenset(v1_candidates[1]) if len(v1_candidates) > 1 else frozenset()
        u1 = frozenset(u1_candidates[0]) if u1_candidates else frozenset()
        u2 = frozenset(u1_candidates[1]) if len(u1_candidates) > 1 else frozenset()
        isolated_edges[(e1, e2)] = (v1, v2, u1, u2)
        new_set = frozenset({(e1, e2)})
        alignment_sets.append(new_set)
        alignment_sets.append(new_set)
        logger.debug(f"Isolated edge processed. Updated alignment sets: {alignment_sets}")
        return alignment_sets, set_matches

    assert f1 is not None and f2 is not None, "[ERROR] Non-isolated edge but no (f1, f2) found!"
    logger.debug(f"Fixed aligned edge (f1, f2): ({f1}, {f2})")
    v1 = frozenset(next((s for s in sets_1 if e1 in s and f1 not in s), set()))
    v2 = frozenset(next((s for s in sets_1 if e1 in s and f1 in s), set()))
    v3 = frozenset(next((s for s in sets_1 if f1 in s and e1 not in s), set()))
    u1 = frozenset(next((s for s in sets_2 if e2 in s and f2 not in s), set()))
    u2 = frozenset(next((s for s in sets_2 if e2 in s and f2 in s), set()))
    u3 = frozenset(next((s for s in sets_2 if f2 in s and e2 not in s), set()))
    logger.debug(f"Unique sets for (f1, f2): v1={v1}, v2={v2}, v3={v3}, u1={u1}, u2={u2}, u3={u3}")
    set_match_v1 = any(v1 in match for match in set_matches.values())
    set_match_v2 = any(v2 in match for match in set_matches.values())
    set_match_v3 = any(v3 in match for match in set_matches.values())
    logger.debug(f"Current set_matches: {set_matches}")
    logger.debug(f"Current alignment_sets: {alignment_sets}")
    if not set_match_v1 and not set_match_v2 and not set_match_v3:
        logger.debug("Case A: No set-match for v1, v2, or v3.")
        singleton_sets = [s for s in alignment_sets if s == frozenset({(f1, f2)})]
        if len(singleton_sets) >= 2:
            w2 = singleton_sets[0]
            w3 = singleton_sets[1]
            logger.debug(f"Found two equal sets for (f1, f2): w2={w2}, w3={w3}")
            updated_w2 = set(w2)
            updated_w2.add((e1, e2))
            alignment_sets[alignment_sets.index(w2)] = frozenset(updated_w2)
            logger.debug(f"Updated w2: {frozenset(updated_w2)} (added (e1, e2))")
            logger.debug(f"w3 remains unchanged: {w3}")
            set_matches[frozenset(updated_w2)] = (v2, u2)
            set_matches[w3] = (v3, u3)
            logger.debug(f"Updated set_matches: w2 -> (v2, u2), w3 -> (v3, u3)")
        else:
            logger.error("Could not find two equal sets for (f1, f2) in alignment_sets.")
        v1_aligned_edges_count = len({e for e in v1 if any(e == edge[0] for edge in current_alignment)})
        logger.debug(f"v1_aligned_edges_count={v1_aligned_edges_count}")
        if v1_aligned_edges_count > 2:
            logger.error("Error: There has to be a matching as there are two edges on e1 in the alignment")
        elif v1_aligned_edges_count == 1:
            logger.debug(f"current_alignment={current_alignment}")
            w1 = frozenset({(e1, e2)})
            alignment_sets.append(w1)
            set_matches[w1] = (v1, u1)
            logger.debug(f"Created new set w1: {w1}")
        elif v1_aligned_edges_count == 2:
            logger.debug("g1g2 Needed!!!!!")
            g1 = next((e for e in v1 if e != e1 and any(e == edge[0] for edge in current_alignment)), None)
            g2 = next((e for e in u1 if e != e2 and any(e == edge[1] for edge in current_alignment)), None)
            if g1 is not None and g2 is not None:
                equal_sets = [s for s in alignment_sets if s == frozenset({(g1, g2)})]
                if len(equal_sets) >= 2:
                    w1 = equal_sets[0]
                    w4 = equal_sets[1]
                    updated_w1 = set(w1)
                    updated_w1.add((e1, e2))
                    alignment_sets[alignment_sets.index(w1)] = frozenset(updated_w1)
                    logger.debug(f"Updated w1: {frozenset(updated_w1)} (added (e1, e2))")
                    logger.debug(f"w4 remains unchanged: {w4}")
                    v0 = frozenset(next((s for s in sets_1 if g1 in s and e1 not in s), set()))
                    u0 = frozenset(next((s for s in sets_2 if g2 in s and e2 not in s), set()))
                    set_matches[frozenset(updated_w1)] = (v1, u1)
                    set_matches[w4] = (v0, u0)
                    logger.debug(f"Updated set_matches: w1 -> (v1, u1), w4 -> (v0, u0)")
        logger.debug(f"Final alignmentAA: {current_alignment}")
        logger.debug(f"Final alignment_sets: {alignment_sets}")
        logger.debug(f"Final set_matches: {set_matches}")
        return alignment_sets, set_matches
    else:
        logger.debug("Case B: Found set-match for v1, v2, or v3.")
        if set_match_v1:
            w1 = next((s3 for s3, (s1, _) in set_matches.items() if s1 == v1), None)
            if w1:
                updated_w1 = set(w1)
                updated_w1.add((e1, e2))
                alignment_sets[alignment_sets.index(w1)] = frozenset(updated_w1)
                set_matches[frozenset(updated_w1)] = set_matches.pop(w1)
                logger.debug(f"Updated w1: {frozenset(updated_w1)} (added (e1, e2))")
            singleton_sets = [s for s in alignment_sets if s == frozenset({(f1, f2)})]
            if len(singleton_sets) >= 2:
                w2 = singleton_sets[0]
                w3 = singleton_sets[1]
                updated_w2 = set(w2)
                updated_w2.add((e1, e2))
                alignment_sets[alignment_sets.index(w2)] = frozenset(updated_w2)
                set_matches[frozenset(updated_w2)] = (v2, u2)
                set_matches[w3] = (v3, u3)
                logger.debug(f"Updated w2: {frozenset(updated_w2)}")
                logger.debug(f"w3 remains unchanged: {w3}")
        elif set_match_v2:
            v1_aligned_edges = {e for e in v1 if any(e == edge[0] for edge in current_alignment)}
            logger.debug(f"v1_aligned_edges(Case B)={v1_aligned_edges}")
            logger.debug(f"current_alignment={current_alignment}")
            w2 = next((s3 for s3, (s1, _) in set_matches.items() if s1 == v2), None)
            if w2:
                updated_w2 = set(w2)
                updated_w2.add((e1, e2))
                alignment_sets[alignment_sets.index(w2)] = frozenset(updated_w2)
                set_matches[frozenset(updated_w2)] = set_matches.pop(w2)
                logger.debug(f"Updated w2: {frozenset(updated_w2)} (added (e1, e2))")
            if len(v1_aligned_edges) == 1:
                w1 = frozenset({(e1, e2)})
                alignment_sets.append(w1)
                set_matches[w1] = (v1, u1)
                logger.debug(f"Created new set w1: {w1}")
            elif len(v1_aligned_edges) == 2:
                g1 = next((e for e in v1 if e != e1 and any(e == edge[0] for edge in current_alignment)), None)
                g2 = next((e for e in u1 if e != e2 and any(e == edge[1] for edge in current_alignment)), None)
                logger.debug(f"We have g1={g1}, g2={g2}")
                if g1 and g2:
                    equal_sets = [s for s in alignment_sets if s == frozenset({(g1, g2)})]
                    logger.debug(f"equal_sets={equal_sets}")
                    if len(equal_sets) >= 2:
                        w1 = equal_sets[0]
                        w4 = equal_sets[1]
                        updated_w1 = set(w1)
                        updated_w1.add((e1, e2))
                        alignment_sets[alignment_sets.index(w1)] = frozenset(updated_w1)
                        set_matches[frozenset(updated_w1)] = (v1, u1)
                        v0 = frozenset(next((s for s in sets_1 if g1 in s and e1 not in s), set()))
                        u0 = frozenset(next((s for s in sets_2 if g2 in s and e2 not in s), set()))
                        set_matches[w4] = (v0, u0)
                        logger.debug(f"Updated w1: {frozenset(updated_w1)}")
                        logger.debug(f"Updated w4: {w4}")
        logger.debug(f"Final alignmentAA: {current_alignment}")
        logger.debug(f"Final alignment_sets: {alignment_sets}")
        logger.debug(f"Final set_matches: {set_matches}")
        return alignment_sets, set_matches


def is_complete_edge(edge: Any) -> bool:
    if isinstance(edge, (tuple, list)):
        return all(is_complete_edge(sub) for sub in edge)
    return edge is not None


def find_best_alignment(
        edges_1: List[Any],
        sets_1: List[Set[Any]],
        edges_2: List[Any],
        sets_2: List[Set[Any]],
        current_depth: int
) -> Tuple[
    List[Tuple[Any, Any]], List[FrozenSet[Tuple[Any, Any]]], int, Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]]]:
    best_alignment: List[Tuple[Any, Any]] = []
    max_cost = 0
    best_alignment_sets: List[FrozenSet[Tuple[Any, Any]]] = []
    best_set_matches: Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]] = {}
    alignment_sets: List[FrozenSet[Tuple[Any, Any]]] = []
    set_matches: Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]] = {}
    isolated_edges: Dict[Tuple[Any, Any], Tuple[Any, Any, Any, Any]] = {}

    def is_feasible_extension(e1: Any, e2: Any, current_alignment: List[Tuple[Any, Any]],
                              sets_1: List[Set[Any]], sets_2: List[Set[Any]]) -> bool:
        if e1 is None or e2 is None:
            return False

        incident_sets_1 = [s for s in sets_1 if e1 in s]
        incident_sets_2 = [s for s in sets_2 if e2 in s]

        if len(incident_sets_1) != 2 or len(incident_sets_2) != 2:
            logger.error("Error: Inputs have to be E-graphs")
            return False

        x, y = incident_sets_1
        a, b = incident_sets_2

        union_neighbors_e1 = set(x).union(set(y)) - {e1}
        union_neighbors_e2 = set(a).union(set(b)) - {e2}

        aligned_neighbors_e1 = {
            (aligned_e1, aligned_e2) for aligned_e1, aligned_e2 in current_alignment if aligned_e1 in union_neighbors_e1
        }
        aligned_neighbors_e2 = {
            (aligned_e1, aligned_e2) for aligned_e1, aligned_e2 in current_alignment if aligned_e2 in union_neighbors_e2
        }
        logger.debug(f"aligned_neighbors_e1: {aligned_neighbors_e1}")
        logger.debug(f"aligned_neighbors_e2: {aligned_neighbors_e2}")
        if aligned_neighbors_e1 != aligned_neighbors_e2:
            logger.debug(f"aligned_neighbors_e1 {aligned_neighbors_e1} != aligned_neighbors_e2 {aligned_neighbors_e2}")
            return False

        x_y_in_set_matches = any(s == x for _, (s, t) in set_matches.items()) and \
                             any(s == y for _, (s, t) in set_matches.items())

        if x_y_in_set_matches:
            logger.debug("Should be already there as an implied edge (for the first coordinate)")
        a_b_in_set_matches = any(t == a for _, (s, t) in set_matches.items()) and \
                             any(t == b for _, (s, t) in set_matches.items())
        if a_b_in_set_matches:
            logger.debug("Should be already there as an implied edge (for the second coordinate)")

        if x_y_in_set_matches or a_b_in_set_matches:
            return False

        return True

    def backtrack_vf2(current_alignment: List[Tuple[Any, Any]],
                      frontier_edges_1: Set[Any],
                      frontier_edges_2: Set[Any],
                      matched_edges_1: Set[Any],
                      matched_edges_2: Set[Any]) -> None:
        nonlocal best_alignment, max_cost, best_set_matches, best_alignment_sets, alignment_sets, set_matches, isolated_edges

        current_cost = sum(1 for edge in current_alignment if is_complete_edge(edge))
        if current_cost > max_cost:
            max_cost = current_cost
            best_alignment = list(current_alignment)
            best_alignment_sets = deepcopy(alignment_sets)
            best_set_matches = deepcopy(set_matches)

        for e1 in list(frontier_edges_1):
            for e2 in list(frontier_edges_2):
                logger.debug(f"EXPLORING candidate edge pair ({e1}, {e2})")
                if is_feasible_extension(e1, e2, current_alignment, sets_1, sets_2):
                    previous_current_alignment = deepcopy(current_alignment)
                    logger.debug(f"Adding edge pair ({e1}, {e2}) to alignment...")
                    current_alignment.append((e1, e2))
                    previous_alignment_sets = deepcopy(alignment_sets)
                    previous_set_matches = deepcopy(set_matches)
                    logger.debug(f"Saved state before adding ({e1}, {e2}):")
                    logger.debug(f"  Alignment Sets: {previous_alignment_sets}")
                    logger.debug(f"  Set Matches: {previous_set_matches}")
                    logger.debug(f"  current_alignment: {previous_current_alignment}")

                    alignment_sets, set_matches = update_alignment_sets(
                        e1, e2, alignment_sets, set_matches, sets_1, sets_2, isolated_edges, current_alignment
                    )

                    set_matches, current_alignment, alignment_sets = check_implied_edge_matches(
                        set_matches, current_alignment, alignment_sets
                    )

                    logger.debug(f"Updated after adding ({e1}, {e2}):")
                    logger.debug(f"  Alignment Sets: {alignment_sets}")
                    logger.debug(f"  Set Matches: {set_matches}")

                    matched_edges_1.add(e1)
                    matched_edges_2.add(e2)
                    frontier_edges_1.remove(e1)
                    frontier_edges_2.remove(e2)

                    new_frontier_edges_1 = {
                                               e for s in sets_1 if any(prev in s for prev in matched_edges_1) for e in
                                               s
                                           } - matched_edges_1
                    new_frontier_edges_2 = {
                                               e for s in sets_2 if any(prev in s for prev in matched_edges_2) for e in
                                               s
                                           } - matched_edges_2

                    backtrack_vf2(
                        current_alignment,
                        frontier_edges_1 | new_frontier_edges_1,
                        frontier_edges_2 | new_frontier_edges_2,
                        matched_edges_1,
                        matched_edges_2,
                    )

                    logger.debug(f"Removing edge pair ({e1}, {e2}) from alignment...")
                    alignment_sets = previous_alignment_sets
                    set_matches = previous_set_matches
                    current_alignment = previous_current_alignment
                    logger.debug(f"Restored state after removing ({e1}, {e2}):")
                    logger.debug(f"  Alignment Sets: {alignment_sets}")
                    logger.debug(f"  Set Matches: {set_matches}")

                    matched_edges_1.remove(e1)
                    matched_edges_2.remove(e2)
                    frontier_edges_1.add(e1)
                    frontier_edges_2.add(e2)
        logger.debug(f"BEST alignment so far: {best_alignment}")

    frontier_edges_1 = set(edges_1)
    frontier_edges_2 = set(edges_2)
    logger.debug(f"Starting backtracking with initial frontier_edges_1: {frontier_edges_1}")
    backtrack_vf2([], frontier_edges_1, frontier_edges_2, set(), set())
    logger.debug(f"BEST alignment after backtracking: {best_alignment}")
    logger.debug(f"BEST set_matches after backtracking: {best_set_matches}")

    best_alignment, best_alignment_sets, best_set_matches = refine_alignment_with_unaligned_edges(
        best_alignment, best_alignment_sets, best_set_matches, sets_1, sets_2, set(edges_1), set(edges_2), current_depth
    )
    logger.debug(f"best_alignment_after refinement: {best_alignment}")
    logger.debug(f"Alignment Sets: {best_alignment_sets}")
    logger.debug(f"Cost: {max_cost}")
    logger.debug(f"Set Matches: {best_set_matches}")

    return best_alignment, best_alignment_sets, max_cost, best_set_matches


def align_multiple_e_graphs(
        edges_list: List[List[Any]],
        sets_list: List[List[Set[Any]]]
) -> Tuple[List[Tuple[Any, Any]], List[FrozenSet[Tuple[Any, Any]]], Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]]]:
    num_graphs = len(edges_list)
    if num_graphs < 2:
        return edges_list[0], sets_list[0], {}

    composite_edges = edges_list[0]
    composite_sets = sets_list[0]
    composite_set_matches: Dict[FrozenSet[Tuple[Any, Any]], Tuple[Any, Any]] = {}
    current_depth = 1
    for i in range(1, num_graphs):
        current_depth += 1
        logger.debug(f"Aligning e-graph number {current_depth}")
        composite_edges, composite_sets, _, composite_set_matches = find_best_alignment(
            composite_edges, composite_sets, edges_list[i], sets_list[i], current_depth
        )
    return composite_edges, composite_sets, composite_set_matches


def load_data_from_json(file_path: str) -> Tuple[List[List[Any]], List[List[Set[Any]]]]:
    """Load edges_list and sets_list from a JSON file and validate the data."""
    with open(file_path, 'r') as f:
        data = json.load(f)

    edges_list = data.get("edges_list", [])
    raw_sets_list = data.get("sets_list", [])

    # Convert each inner list in sets_list to a list of Python sets.
    sets_list: List[List[Set[Any]]] = []
    for graph in raw_sets_list:
        # Each graph is a list of collections; convert each to a set.
        sets_list.append([set(item) for item in graph])

    # Check for distinct edge labels across e-graphs.
    global_edge_to_graph: Dict[Any, int] = {}
    for i, graph_edges in enumerate(edges_list):
        for edge in graph_edges:
            if edge in global_edge_to_graph:
                raise ValueError(
                    f"Edge label '{edge}' appears in both e-graph {global_edge_to_graph[edge]} and e-graph {i}. "
                    "It is recommended to use distinct edge labels for clearer output and interpretation."
                )
            global_edge_to_graph[edge] = i

    # For each e-graph, ensure each edge appears in exactly 0 or 2 sets.
    for i, (graph_edges, graph_sets) in enumerate(zip(edges_list, sets_list)):
        for edge in graph_edges:
            appearance_count = sum(1 for s in graph_sets if edge in s)
            if appearance_count == 0:
                # Automatically add two singleton sets {edge} and {edge}.
                logger.debug(f"In e-graph {i}, edge '{edge}' appears 0 times. Adding two singleton sets for it.")
                graph_sets.append({edge})
                graph_sets.append({edge})
            elif appearance_count != 2:
                raise ValueError(
                    f"In e-graph {i}, edge '{edge}' appears {appearance_count} times. "
                    "Each edge must appear in either 0 or 2 sets. Please fix your input."
                )
    return edges_list, sets_list


def main() -> None:
    parser = argparse.ArgumentParser(description="Align multiple E-graphs.")
    parser.add_argument(
        "--log-level", type=str, default="INFO",
        help="Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL). Default is INFO."
    )
    parser.add_argument(
        "--log-file", type=str, default=None,
        help="Optional file to output logs."
    )
    parser.add_argument(
        "--input", type=str, default=None,
        help="Path to JSON file containing 'edges_list' and 'sets_list'."
    )
    args = parser.parse_args()
    configure_logging(args.log_level, args.log_file)

    # Load data either from a JSON file or use fallback defaults.
    if args.input:
        try:
            edges_list, sets_list = load_data_from_json(args.input)
            logger.info(f"Loaded input from {args.input}")
        except Exception as e:
            logger.error(f"Failed to load input data: {e}")
            return
    else:
        # Fallback default inputs.
        edges_list = [
            ['e1', 'e2', 'e3', 'e4'],
            ['g1', 'g2', 'g3', 'g4'],
            ['f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7'],
            ['k1', 'k2', 'k3', 'k4']
        ]
        sets_list = [
            [{'e1', 'e2'}, {'e2', 'e3'}, {'e4', 'e3'}, {'e1', 'e4'}],
            [{'g1'}, {'g2', 'g3'}, {'g4', 'g3'}, {'g1', 'g2'}, {'g4'}],
            [{'f1', 'f2', 'f5'}, {'f2', 'f3'}, {'f3', 'f4'}, {'f4', 'f1'}, {'f5', 'f6'}, {'f6', 'f7'}, {'f7'}],
            [{'k1'}, {'k2', 'k3'}, {'k4', 'k3'}, {'k1', 'k2'}, {'k4'}]
        ]
        logger.info("No input file provided. Using default test data.")

    # Proceed with alignment.
    alignment_edges, alignment_sets, set_matches = align_multiple_e_graphs(edges_list, sets_list)

    logger.info("Final Alignment Edges:")
    for edge in alignment_edges:
        logger.info(edge)

    logger.info("Final Alignment Sets:")
    for idx, edge_set in enumerate(alignment_sets):
        logger.info(f"Set {idx + 1}: {edge_set}")

    logger.info("Final Set Matches:")
    for s3, (s1, s2) in set_matches.items():
        logger.info(f"{s3} -> {s1}, {s2}")


if __name__ == "__main__":
    main()
