import argparse
import logging
import json
from collections import defaultdict, deque, Counter
from itertools import combinations
from copy import deepcopy
from typing import Any, List, Set, Tuple, Dict, Optional, FrozenSet, Union

import networkx as nx
import matplotlib.pyplot as plt

CompositeEdge = Union[Any, Tuple[Any, ...]]


# ----------------- Logging Configuration -----------------
def configure_logging(log_level: str, log_file: Optional[str] = None) -> None:
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
    global logger
    logger = logging.getLogger(__name__)


# ----------------- Padding Helpers -----------------
def pad_none(depth: int) -> Optional[Any]:
    """
    Recursively create a padded structure based on the specified depth.
    This helper is used when no edge is present.
    """
    if depth <= 1:
        return None
    else:
        return (pad_none(depth - 1), None)


def pad_edge(e: CompositeEdge, current_depth: int) -> CompositeEdge:
    """
    Wraps an edge 'e' in a padded structure according to the current depth.
    When one side of an alignment is missing, pad_edge creates a structure
    matching the complete edge format.
    """
    if current_depth <= 1:
        return e
    else:
        return (pad_none(current_depth - 1), e)


# ----------------- Alignment Utilities -----------------
def is_complete_edge(edge: CompositeEdge) -> bool:
    """
    Recursively checks if the given edge (or nested structure of edges) is complete,
    i.e. that no part of the edge is None.
    """
    if isinstance(edge, (tuple, list)):
        return all(is_complete_edge(sub_edge) for sub_edge in edge)
    return edge is not None

SetMatchValue = Union[
    Tuple[FrozenSet[Any], FrozenSet[Any]],
    Tuple[str, FrozenSet[Any]],
    Tuple[FrozenSet[Any], str]
]

def check_implied_edge_matches(
        set_matches: Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],
        alignment_edges: List[CompositeEdge],
        alignment_sets: List[FrozenSet[CompositeEdge]],
        edge_type_map_1: Dict[CompositeEdge, str],
        edge_type_map_2: Dict[CompositeEdge, str],
        set_type_map_1: Dict[FrozenSet[Any], str],
        set_type_map_2: Dict[FrozenSet[Any], str],
        alignment_edge_type_map: Optional[Dict[str, str]] = None,
        alignment_set_type_map: Optional[Dict[FrozenSet[CompositeEdge], str]] = None
) -> Tuple[
    Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],
    List[CompositeEdge],
    List[FrozenSet[CompositeEdge]],
    Optional[Dict[str, str]],
    Optional[Dict[FrozenSet[CompositeEdge], str]]
]:
    logger.debug("Checking for implied edge matches with type constraints on edges and sets...")
    updated = True
    while updated:
        updated = False
        for (w1_entry, (v1, u1)), (w2_entry, (v2, u2)) in combinations(set_matches.items(), 2):
            w1 = w1_entry[0]
            w2 = w2_entry[0]
            if w1 != w2:
                # Ensure the corresponding sets have matching types.
                type_v1 = set_type_map_1.get(v1)
                type_v2 = set_type_map_1.get(v2)
                type_u1 = set_type_map_2.get(u1)
                type_u2 = set_type_map_2.get(u2)
                logger.debug(
                    f"type_v1:{type_v1}== type_u1:{type_u1}, type_v2:{type_v2}== type_u2:{type_u2}")

                if type_v1 != type_u1 or type_v2 != type_u2:
                    logger.debug(
                        f"Skipping pair {w1} and {w2} due to mismatched set types: {type_v1} vs {type_u1} or {type_v2} vs {type_u2}")
                    continue
                v1_set, u1_set = frozenset(v1), frozenset(u1)
                v2_set, u2_set = frozenset(v2), frozenset(u2)
                v_intersection = v1_set.intersection(v2_set)
                u_intersection = u1_set.intersection(u2_set)
                logger.debug(f"v_intersection={v_intersection}")
                logger.debug(f"u_intersection={u_intersection}")
                if v_intersection and u_intersection:
                    e1 = e2 = None
                    found_candidate = False
                    for candidate_e1 in v_intersection:
                        for candidate_e2 in u_intersection:
                            logger.debug(f"candidate_e1={candidate_e1}")
                            logger.debug(f"candidate_e2={candidate_e2}")
                            type1 = edge_type_map_1.get(str(candidate_e1))
                            type2 = edge_type_map_2.get(str(candidate_e2))
                            logger.debug(f"type1={type1}")
                            logger.debug(f"type2={type2}")
                            # If both types are known but do not match, this is a forbidden implied edge.
                            if type1 is not None and type2 is not None and type1 != type2:
                                logger.debug(
                                    f"Forbidden implied edge: candidate pair ({candidate_e1}, {candidate_e2}) has mismatched types ({type1} vs {type2}). Pruning branch.")
                                raise Exception("Forbidden implied edge encountered: alignment branch pruned.")
                            # If types match, record the candidate pair.
                            if type1 is not None and type1 == type2:
                                e1 = candidate_e1
                                e2 = candidate_e2
                                found_candidate = True
                                break
                        if found_candidate:
                            break
                    if found_candidate:
                        logger.debug(f"alignment_edges: {alignment_edges}")
                        if (e1, e2) not in alignment_edges:
                        #if not any((e1, e2) in s for s in alignment_sets):
                            logger.debug(f"Found implied edge: ({e1}, {e2})")
                            logger.debug(f"Found implied edge: ({e1}, {e2}) with edge type '{edge_type_map_1[str(e1)]}'")
                            alignment_edges.append((e1, e2))
                            alignment_edge_type_map[str((e1, e2))] = edge_type_map_1[str(e1)]
                            updated_w1 = set(w1)
                            updated_w1.add((e1, e2))
                            new_w1 = frozenset(updated_w1)
                            updated_w2 = set(w2)
                            updated_w2.add((e1, e2))
                            new_w2 = frozenset(updated_w2)
                            alignment_sets[alignment_sets.index(w1)] = new_w1
                            alignment_sets[alignment_sets.index(w2)] = new_w2
                            alignment_set_type_map[new_w1] = type_v1
                            alignment_set_type_map[new_w2] = type_v2

                            if (w1, type_v1) in set_matches:
                                set_matches[(new_w1, type_v1)] = set_matches.pop((w1, type_v1))
                            if (w2, type_v2) in set_matches:
                                set_matches[(new_w2, type_v2)] = set_matches.pop((w2, type_v2))
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
    return set_matches, alignment_edges, alignment_sets, alignment_edge_type_map, alignment_set_type_map


def update_alignment_sets(
        e1: CompositeEdge,
        e2: CompositeEdge,
        alignment_sets: List[FrozenSet[CompositeEdge]],
        set_matches: Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],
        sets_1: List[Set[Any]],
        sets_2: List[Set[Any]],
        current_alignment: List[CompositeEdge],
        edge_type_map_1: Dict[CompositeEdge, str],
        edge_type_map_2: Dict[CompositeEdge, str],
        set_type_map_1: Dict[FrozenSet[Any], str],
        set_type_map_2: Dict[FrozenSet[Any], str],
        alignment_edge_type_map: Dict[str, str],
        alignment_set_type_map: Optional[Dict[FrozenSet[CompositeEdge], str]] = None
) -> Tuple[
    List[FrozenSet[CompositeEdge]],
    Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],
    Dict[str, str],
    Optional[Dict[FrozenSet[CompositeEdge], str]]
]:
    # Ensure alignment_set_type_map is initialized.
    if alignment_set_type_map is None:
        alignment_set_type_map = {}

    logger.debug(f"Updating alignment sets for edge ({e1}, {e2}) with type constraints...")
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
        logger.debug("Edge is isolated.")
        # Expect exactly two candidate sets from each graph.
        if len(sets_1) == 0 or len(sets_2) == 0:
            logger.critical("Insufficient candidate sets: no candidate found.")
            raise ValueError("Insufficient candidate sets: no candidate found for one of the graphs.")

        v1_candidates = [s for s in sets_1 if e1 in s]
        u1_candidates = [s for s in sets_2 if e2 in s]
        if len(v1_candidates) < 2 or len(u1_candidates) < 2:
            logger.critical("Insufficient candidate sets: fewer than 2 candidates found.")
            raise ValueError("Insufficient candidate sets: each edge must appear in exactly 2 sets.")
        v1 = frozenset(v1_candidates[0])
        v2 = frozenset(v1_candidates[1])
        u1 = frozenset(u1_candidates[0])
        u2 = frozenset(u1_candidates[1])
        # Retrieve types from the provided maps.
        type_v1 = set_type_map_1.get(v1)
        type_v2 = set_type_map_1.get(v2)
        type_u1 = set_type_map_2.get(u1)
        type_u2 = set_type_map_2.get(u2)
        logger.debug(f"Candidate types: v1={type_v1}, v2={type_v2}, u1={type_u1}, u2={type_u2}")

        ## If the candidate sets in each graph do not differ, we fall back.
        if type_v1 == type_v2 and type_u1 == type_u2:
            forward = False
            backward = False

            new_set1 = frozenset({(e1, e2)})
            new_set2 = frozenset({(e1, e2)})
            label1='--An_Implicit_Label--'
            label2='--Another_Implicit_Label--'
            alignment_set_type_map[new_set1] = label1
            alignment_set_type_map[new_set2] =  label2
            # Forward direction:
            set_matches[(new_set1, label1)] = (v1, u1)
            set_matches[(new_set2, label2)] = (v2, u2)

            # now we check if it creates any forbidden implied edge
            for (w1_entry, (y1, z1)), (w2_entry, (y2, z2)) in combinations(set_matches.items(), 2):
                w1 = w1_entry[0]
                w2 = w2_entry[0]
                if w1 != w2:
                    # Ensure the corresponding sets have matching types.
                    type_y1 = set_type_map_1.get(y1)
                    type_y2 = set_type_map_1.get(y2)
                    type_z1 = set_type_map_2.get(z1)
                    type_z2 = set_type_map_2.get(z2)

                    if type_y1 != type_z1 or type_y2 != type_z2:
                        logger.debug(
                            f"Skipping pair {w1} and {w2} due to mismatched set types: {type_v1} vs {type_u1} or {type_v2} vs {type_u2}")
                        continue
                    y1_set, z1_set = frozenset(y1), frozenset(z1)
                    y2_set, z2_set = frozenset(y2), frozenset(z2)
                    y_intersection = y1_set.intersection(y2_set)
                    z_intersection = z1_set.intersection(z2_set)
                    if y_intersection and z_intersection:
                        for candidate_e1 in y_intersection:
                            for candidate_e2 in z_intersection:
                                type1 = edge_type_map_1.get(str(candidate_e1))
                                type2 = edge_type_map_2.get(str(candidate_e2))
                                # If both types are known but do not match, this direction (forward) yields a forbidden implied edge.
                                if type1 is not None and type2 is not None and type1 != type2:
                                    forward = True
                                    logger.debug(
                                        f"forward set_matches for the isolated edge will result in ({candidate_e1}, {candidate_e2}) with mismatched types ({type1} vs {type2}) which is a forbidden implied edge.")

            # Backward direction: First we remove
            del set_matches[(new_set1, label1)]
            del set_matches[(new_set2, label2)]
            labela = '--An_Implicit_Label_A--'
            labelb = '--An_Implicit_Label_B--'
            set_matches[(new_set1, labela)] = (v1, u2)
            set_matches[(new_set2, labelb)] = (v2, u1)
            # now we check if the backward direction creates any forbidden implied edge
            for (w1_entry, (y1, z1)), (w2_entry, (y2, z2)) in combinations(set_matches.items(), 2):
                w1 = w1_entry[0]
                w2 = w2_entry[0]
                if w1 != w2:
                    # Ensure the corresponding sets have matching types.
                    type_y1 = set_type_map_1.get(y1)
                    type_y2 = set_type_map_1.get(y2)
                    type_z1 = set_type_map_2.get(z1)
                    type_z2 = set_type_map_2.get(z2)

                    if type_y1 != type_z1 or type_y2 != type_z2:
                        logger.debug(
                            f"Skipping pair {w1} and {w2} due to mismatched set types: {type_v1} vs {type_u1} or {type_v2} vs {type_u2}")
                        continue
                    y1_set, z1_set = frozenset(y1), frozenset(z1)
                    y2_set, z2_set = frozenset(y2), frozenset(z2)
                    y_intersection = y1_set.intersection(y2_set)
                    z_intersection = z1_set.intersection(z2_set)
                    if y_intersection and z_intersection:
                        for candidate_e1 in y_intersection:
                            for candidate_e2 in z_intersection:
                                type1 = edge_type_map_1.get(str(candidate_e1))
                                type2 = edge_type_map_2.get(str(candidate_e2))
                                # If both types are known but do not match, this direction (forward) yields a forbidden implied edge.
                                if type1 is not None and type2 is not None and type1 != type2:
                                    backward = True
                                    logger.debug(
                                        f"Backward set_matches for the isolated edge will result in ({candidate_e1}, {candidate_e2}) with mismatched types ({type1} vs {type2}) which is a forbidden implied edge.")
            del set_matches[(new_set1, labela)]
            del set_matches[(new_set2, labelb)]

            if forward and backward:
                logger.debug(
                    f"Any set_matches for the isolated edge ({e1}, {e2}) will result in a forbidden implied edege: alignment branch pruned.")
                raise Exception("Forbidden implied edge encountered: alignment branch pruned.")

            logger.debug(
                f"Isolated edge {e1}/{e2}: Candidate sets do not differ in type. Deferring decision for Set_matches.")
            # Fall back to default handling.
            new_set = frozenset({(e1, e2)})
            alignment_sets.append(new_set)
            alignment_sets.append(new_set)
            alignment_set_type_map[new_set] = set_type_map_1.get(frozenset(v1_candidates[0]), type_v1)
            return alignment_sets, set_matches, alignment_edge_type_map, alignment_set_type_map

        # Now check for valid pairing possibilities.
        if type_v1 == type_u1 and type_v2 == type_u2:
            # Possibility 1: v1 pairs with u1, v2 with u2.
            new_set1 = frozenset({(e1, e2)})
            new_set2 = frozenset({(e1, e2)})
            alignment_sets.append(new_set1)
            alignment_sets.append(new_set2)
            # Use composite keys: (alignment_set, type)
            set_matches[(new_set1, type_v1)] = (v1, u1)
            set_matches[(new_set2, type_v2)] = (v2, u2)
            alignment_set_type_map[new_set1] = type_v1
            alignment_set_type_map[new_set2] = type_v2
            logger.debug(f"Created two alignment sets: {new_set1} -> {type_v1}, {new_set2} -> {type_v2}")
            return alignment_sets, set_matches, alignment_edge_type_map, alignment_set_type_map
        elif type_v1 == type_u2 and type_v2 == type_u1:
            # Possibility 2: v1 pairs with u2, v2 with u1.
            new_set1 = frozenset({(e1, e2)})
            new_set2 = frozenset({(e1, e2)})
            alignment_sets.append(new_set1)
            alignment_sets.append(new_set2)
            set_matches[(new_set1, type_v1)] = (v1, u2)
            set_matches[(new_set2, type_v2)] = (v2, u1)
            alignment_set_type_map[new_set1] = type_v1
            alignment_set_type_map[new_set2] = type_v2
            logger.debug(f"Created two alignment sets: {new_set1} -> {type_v1}, {new_set2} -> {type_v2}")
            return alignment_sets, set_matches, alignment_edge_type_map, alignment_set_type_map
        else:
            logger.critical("Unexpected candidate endpoint type configuration encountered. This should not happen!")
            raise ValueError("Unexpected candidate endpoint type configuration.")

    assert f1 is not None and f2 is not None, "[ERROR] Non-isolated edge but no (f1, f2) found!"
    logger.debug(f"Fixed aligned edge (f1, f2): ({f1}, {f2})")
    v1 = frozenset(next((s for s in sets_1 if e1 in s and f1 not in s), set()))
    v2 = frozenset(next((s for s in sets_1 if e1 in s and f1 in s), set()))
    v3 = frozenset(next((s for s in sets_1 if f1 in s and e1 not in s), set()))
    u1 = frozenset(next((s for s in sets_2 if e2 in s and f2 not in s), set()))
    u2 = frozenset(next((s for s in sets_2 if e2 in s and f2 in s), set()))
    u3 = frozenset(next((s for s in sets_2 if f2 in s and e2 not in s), set()))
    type_v1 = set_type_map_1.get(v1)
    type_v2 = set_type_map_1.get(v2)
    type_v3 = set_type_map_1.get(v3)
    
    if set_type_map_1.get(v1) != set_type_map_2.get(u1) or set_type_map_1.get(v2) != set_type_map_2.get(u2) \
            or set_type_map_1.get(v3) != set_type_map_2.get(u3):
        logger.debug(
            f"Candidate sets have mismatched types: {set_type_map_1.get(v1)}, {set_type_map_1.get(v2)}, {set_type_map_1.get(v3)}")
        logger.debug(
            f"Candidate sets have mismatched types: {set_type_map_2.get(u1)}, {set_type_map_2.get(u2)}, {set_type_map_1.get(u3)}")
        logger.debug(
            f"Adding ({e1}, {e2}) is inconsistent with the current set_match")
        raise Exception("Inconsistent with the set_match. alignment branch pruned.")

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
            new_w2 = frozenset(updated_w2)
            alignment_sets[alignment_sets.index(w2)] = new_w2
            alignment_set_type_map[new_w2] = type_v2
            logger.debug(f"Updated w2: {frozenset(updated_w2)} (added (e1, e2))")
            logger.debug(f"w3 remains unchanged: {w3}")
            set_matches[(new_w2, type_v2)] = (v2, u2)
            set_matches[(w3, type_v3)] = (v3, u3)
            logger.debug(f"Updated set_matches: w2 -> (v2, u2), w3 -> (v3, u3)")
        else:
            logger.error("Could not find two equal sets for (f1, f2) in alignment_sets.")
        v1_aligned_edges_count = len({e for e in v1 if any(e == edge[0] for edge in current_alignment)})
        logger.debug(f"v1_aligned_edges_count={v1_aligned_edges_count}")
        if v1_aligned_edges_count > 2:
            logger.error("Error: There are more than two aligned edges on v1 in the alignment")
        elif v1_aligned_edges_count == 1:
            logger.debug(f"current_alignment={current_alignment}")
            w1 = frozenset({(e1, e2)})
            alignment_sets.append(w1)
            alignment_set_type_map[w1] = type_v1
            set_matches[(w1,type_v1)] = (v1, u1)
            logger.debug(f"Created new set w1: {w1}")
        elif v1_aligned_edges_count == 2:
            logger.debug("Additional matching required (g1g2 case)!")
            g1 = next((e for e in v1 if e != e1 and any(e == edge[0] for edge in current_alignment)), None)
            g2 = next((e for e in u1 if e != e2 and any(e == edge[1] for edge in current_alignment)), None)
            if g1 is not None and g2 is not None:
                equal_sets = [s for s in alignment_sets if s == frozenset({(g1, g2)})]
                if len(equal_sets) >= 2:
                    w1 = equal_sets[0]
                    w4 = equal_sets[1]
                    updated_w1 = set(w1)
                    updated_w1.add((e1, e2))
                    new_w1 = frozenset(updated_w1)
                    alignment_sets[alignment_sets.index(w1)] = new_w1
                    alignment_set_type_map[new_w1] = type_v1
                    logger.debug(f"Updated w1: {frozenset(updated_w1)} (added (e1, e2))")
                    logger.debug(f"w4 remains unchanged: {w4}")
                    v0 = frozenset(next((s for s in sets_1 if g1 in s and e1 not in s), set()))
                    u0 = frozenset(next((s for s in sets_2 if g2 in s and e2 not in s), set()))
                    type_v0 = set_type_map_1.get(v0)
                    set_matches[(new_w1, type_v1)] = (v1, u1)
                    set_matches[(w4, type_v0)] = (v0, u0)
                    logger.debug(f"Updated set_matches: w1 -> (v1, u1), w4 -> (v0, u0)")
    else:
        logger.debug("Case B: Found existing set-match for v1, v2, or v3.")
        if set_match_v1:
            w1_entry = next((s3 for s3, (s1, _) in set_matches.items() if s1 == v1), None)
            if w1_entry:
                w1 = w1_entry[0]
                updated_w1 = set(w1)
                updated_w1.add((e1, e2))
                new_w1 = frozenset(updated_w1)
                alignment_sets[alignment_sets.index(w1)] = new_w1
                alignment_set_type_map[new_w1] = type_v1
                set_matches[(new_w1, type_v1)] = set_matches.pop((w1, type_v1))
                logger.debug(f"Updated w1: {frozenset(updated_w1)} (added (e1, e2))")
            singleton_sets = [s for s in alignment_sets if s == frozenset({(f1, f2)})]
            if len(singleton_sets) >= 2:
                w2 = singleton_sets[0]
                w3 = singleton_sets[1]
                updated_w2 = set(w2)
                updated_w2.add((e1, e2))
                new_w2 = frozenset(updated_w2)
                alignment_sets[alignment_sets.index(w2)] = new_w2
                alignment_set_type_map[new_w2] = type_v2
                set_matches[(new_w2, type_v2)] = (v2, u2)
                set_matches[(w3, type_v3)] = (v3, u3)
                logger.debug(f"Updated w2: {frozenset(updated_w2)}")
                logger.debug(f"w3 remains unchanged: {w3}")
        elif set_match_v2:
            v1_aligned_edges = {e for e in v1 if any(e == edge[0] for edge in current_alignment)}
            logger.debug(f"v1_aligned_edges (Case B)={v1_aligned_edges}")
            logger.debug(f"current_alignment={current_alignment}")
            w2_entry = next((s3 for s3, (s1, _) in set_matches.items() if s1 == v2), None)
            if w2_entry:
                w2 = w2_entry[0]
                updated_w2 = set(w2)
                logger.debug(f"w2={updated_w2}")
                updated_w2.add((e1, e2))
                new_w2 = frozenset(updated_w2)
                alignment_sets[alignment_sets.index(w2)] = new_w2
                alignment_set_type_map[new_w2] = type_v2
                set_matches[(new_w2, type_v2)] = set_matches.pop((w2,type_v2))
                logger.debug(f"Updated w2: {frozenset(updated_w2)} (added (e1, e2))")
            if len(v1_aligned_edges) == 1:
                w1 = frozenset({(e1, e2)})
                alignment_sets.append(w1)
                alignment_set_type_map[w1] = type_v1
                set_matches[(w1, type_v1)] = (v1, u1)
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
                        new_w1 = frozenset(updated_w1)
                        alignment_sets[alignment_sets.index(w1)] = new_w1
                        alignment_set_type_map[new_w1] = type_v1
                        set_matches[(new_w1, type_v1)] = (v1, u1)
                        v0 = frozenset(next((s for s in sets_1 if g1 in s and e1 not in s), set()))
                        u0 = frozenset(next((s for s in sets_2 if g2 in s and e2 not in s), set()))
                        type_v0 = set_type_map_1.get(v0)
                        set_matches[(w4, type_v0)] = (v0, u0)
                        logger.debug(f"Updated set_matches: w1 -> (v1, u1), w4 -> (v0, u0)")
    logger.debug(f"Final alignment (current_alignment): {current_alignment}")
    logger.debug(f"Final alignment_sets: {alignment_sets}")
    logger.debug(f"Final set_matches: {set_matches}")
    return alignment_sets, set_matches, alignment_edge_type_map, alignment_set_type_map


def handling_set_matches_isolated_edges(
        current_alignment: List[CompositeEdge],
        alignment_sets: List[FrozenSet[CompositeEdge]],
        set_matches: Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],
        alignment_edge_type_map: Dict[str, str],
        alignment_set_type_map: Dict[FrozenSet[CompositeEdge], str],
        sets_1: List[Set[Any]],
        sets_2: List[Set[Any]],
        edge_type_map_1: Dict[CompositeEdge, str],
        edge_type_map_2: Dict[CompositeEdge, str],
        set_type_map_1: Dict[FrozenSet[Any], str],
        set_type_map_2: Dict[FrozenSet[Any], str]
) -> Tuple[
    List[CompositeEdge],
    List[FrozenSet[CompositeEdge]],
    Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],
    Dict[str, str],
    Dict[FrozenSet[CompositeEdge], str]
]:


    logger.debug("Handling isolated edges in set_matches (undecided direction).")

    def edge_already_decided(edge, set_matches):
        """
        Returns True if there is any key in set_matches, whose first coordinate is a singleton set,
        where the only element has its first component equal to the first element of edge
        or its second component equal to the second element of edge.
        """
        e1, e2 = edge
        for (alignment_set, _), (v_candidate, u_candidate) in set_matches.items():
            if len(alignment_set) == 1:
                candidate = next(iter(alignment_set))
                # Here we check if the candidate's first coordinate equals e1
                # or candidate's second coordinate equals e2.
                if candidate[0] == e1 or candidate[1] == e2:
                    return True
        return False
    # First, extract candidate isolated edges.
    # We consider an alignment set isolated if it is a singleton.
    isolated_edge_sets = [s for s in alignment_sets if len(s) == 1]
    # Count frequencies: we need those that appear exactly twice.
    edge_freq = {}
    for s in isolated_edge_sets:
        edge = next(iter(s))
        edge_freq[edge] = edge_freq.get(edge, 0) + 1
    # I will be the set of isolated edges appearing exactly twice.
    I = {edge for edge, count in edge_freq.items() if count == 2 and not edge_already_decided(edge, set_matches)}
    logger.debug(f"Isolated edges extracted for handling (after filtering): {I}")

    # For each isolated edge, try to determine a safe direction.
    for edge in I:
        # Assume edge is a tuple (a, b)
        a, b = edge
        # Retrieve candidate endpoint sets from the respective graphs.
        v1_candidates = [s for s in sets_1 if a in s]
        u1_candidates = [s for s in sets_2 if b in s]
        if len(v1_candidates) < 2 or len(u1_candidates) < 2:
            logger.critical(f"Insufficient candidate sets for edge {edge}.")
            raise ValueError("Insufficient candidate sets: each edge must appear in exactly 2 sets.")

        # For simplicity, pick the first two candidates for each graph.
        v1 = frozenset(v1_candidates[0])
        v2 = frozenset(v1_candidates[1])
        u1 = frozenset(u1_candidates[0])
        u2 = frozenset(u1_candidates[1])

        type_v1 = set_type_map_1.get(v1)
        type_v2 = set_type_map_1.get(v2)
        type_u1 = set_type_map_2.get(u1)
        type_u2 = set_type_map_2.get(u2)

        logger.debug(
            f"For isolated edge {edge}: Candidate types: v1={type_v1}, v2={type_v2}, u1={type_u1}, u2={type_u2}")

        # We now want to test two possible pairing directions.
        forward_forbidden = False
        backward_forbidden = False

        # Prepare two new singleton sets (they will be used as keys in set_matches).
        new_set_forward = frozenset({edge})
        new_set_backward = frozenset({edge})

        # -- Forward Direction: Pair v1 with u1 and v2 with u2.
        # Temporarily add these entries.
        set_matches[(new_set_forward, 'forward')] = (v1, u1)
        set_matches[(new_set_forward, 'forward2')] = (v2, u2)
        # Check all pairings in set_matches for forbidden implied edges.
        for (w1_entry, (y1, z1)), (w2_entry, (y2, z2)) in combinations(set_matches.items(), 2):
            if w1_entry[0] != w2_entry[0]:
                type_y1 = set_type_map_1.get(y1)
                type_y2 = set_type_map_1.get(y2)
                type_z1 = set_type_map_2.get(z1)
                type_z2 = set_type_map_2.get(z2)
                # Skip if types do not match already.
                if type_y1 != type_z1 or type_y2 != type_z2:
                    continue
                y1_set, z1_set = frozenset(y1), frozenset(z1)
                y2_set, z2_set = frozenset(y2), frozenset(z2)
                if y1_set.intersection(y2_set) and z1_set.intersection(z2_set):
                    # For each candidate edge in the intersection, check types.
                    for candidate_e1 in y1_set.intersection(y2_set):
                        for candidate_e2 in z1_set.intersection(z2_set):
                            type1 = edge_type_map_1.get(str(candidate_e1))
                            type2 = edge_type_map_2.get(str(candidate_e2))
                            if type1 is not None and type2 is not None and type1 != type2:
                                forward_forbidden = True
                                logger.debug(
                                    f"Forward direction for isolated edge {edge} forbidden: candidate ({candidate_e1}, {candidate_e2}) mismatched ({type1} vs {type2}).")
        # Clean up forward temporary entries.
        del set_matches[(new_set_forward, 'forward')]
        del set_matches[(new_set_forward, 'forward2')]

        # -- Backward Direction: Pair v1 with u2 and v2 with u1.
        set_matches[(new_set_backward, 'backward')] = (v1, u2)
        set_matches[(new_set_backward, 'backward2')] = (v2, u1)
        for (w1_entry, (y1, z1)), (w2_entry, (y2, z2)) in combinations(set_matches.items(), 2):
            if w1_entry[0] != w2_entry[0]:
                type_y1 = set_type_map_1.get(y1)
                type_y2 = set_type_map_1.get(y2)
                type_z1 = set_type_map_2.get(z1)
                type_z2 = set_type_map_2.get(z2)
                if type_y1 != type_z1 or type_y2 != type_z2:
                    continue
                y1_set, z1_set = frozenset(y1), frozenset(z1)
                y2_set, z2_set = frozenset(y2), frozenset(z2)
                if y1_set.intersection(y2_set) and z1_set.intersection(z2_set):
                    for candidate_e1 in y1_set.intersection(y2_set):
                        for candidate_e2 in z1_set.intersection(z2_set):
                            type1 = edge_type_map_1.get(str(candidate_e1))
                            type2 = edge_type_map_2.get(str(candidate_e2))
                            if type1 is not None and type2 is not None and type1 != type2:
                                backward_forbidden = True
                                logger.debug(
                                    f"Backward direction for isolated edge {edge} forbidden: candidate ({candidate_e1}, {candidate_e2}) mismatched ({type1} vs {type2}).")
        del set_matches[(new_set_backward, 'backward')]
        del set_matches[(new_set_backward, 'backward2')]

        # If both directions are forbidden, we prune the branch.
        if forward_forbidden and backward_forbidden:
            logger.debug(f"Both directions for isolated edge {edge} yield forbidden implied edges. Pruning branch.")
            raise Exception(
                "Forbidden implied edge encountered in handling isolated set_matches: alignment branch pruned.")

        # Otherwise, select the safe direction.
        # Here, we choose forward if safe; otherwise backward.
        if not forward_forbidden:
            safe_direction = 'forward'
            set_matches[(new_set_forward, 'forward')] = (v1, u1)
            set_matches[(new_set_forward, 'forward2')] = (v2, u2)
            # Also, add the corresponding new singleton set into alignment_sets if not already present.
            if new_set_forward not in alignment_sets:
                alignment_sets.append(new_set_forward)
            logger.debug(f"Isolated edge {edge} resolved using forward direction.")
        elif not backward_forbidden:
            safe_direction = 'backward'
            set_matches[(new_set_backward, 'backward')] = (v1, u2)
            set_matches[(new_set_backward, 'backward2')] = (v2, u1)
            if new_set_backward not in alignment_sets:
                alignment_sets.append(new_set_backward)
            logger.debug(f"Isolated edge {edge} resolved using backward direction.")
        else:
            # This branch is unreachable because of the earlier check.
            logger.critical("Unexpected state: no safe direction available despite check.")
            raise Exception("Unexpected error in handling isolated edges.")

        # After deciding on a safe direction for this isolated edge, we call check_implied_edge_matches
        # to propagate any additional consequences.
        try:
            set_matches, current_alignment, alignment_sets, alignment_edge_type_map, alignment_set_type_map = check_implied_edge_matches(
                set_matches, current_alignment, alignment_sets,
                edge_type_map_1, edge_type_map_2, set_type_map_1, set_type_map_2,
                alignment_edge_type_map, alignment_set_type_map
            )
        except Exception as ex:
            logger.debug(f"Forbidden implied edge detected after handling isolated edge {edge}: {ex}. Pruning branch.")
            raise

    logger.debug("Completed handling of isolated set_matches edges.")
    return current_alignment, alignment_sets, set_matches, alignment_edge_type_map, alignment_set_type_map


def find_best_alignment(
        edges_1: List[CompositeEdge],
        sets_1: List[Set[Any]],
        edges_2: List[CompositeEdge],
        sets_2: List[Set[Any]],
        current_depth: int,
        edge_type_map_1: Dict[CompositeEdge, str],
        edge_type_map_2: Dict[CompositeEdge, str],
        set_type_map_1: Dict[FrozenSet[Any], str],
        set_type_map_2: Dict[FrozenSet[Any], str]
) -> Tuple[
    List[CompositeEdge],                                   # alignment_edges
    List[FrozenSet[CompositeEdge]],                        # alignment_sets
    int,                                                     # max_cost
    Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],   # set_matches
    Dict[str, str],                              # alignment_edge_type_map
    Dict[FrozenSet[CompositeEdge], str]                     # alignment_set_type_map
]:
    best_alignment: List[CompositeEdge] = []
    best_alignment_sets: List[FrozenSet[CompositeEdge]] = []
    best_set_matches: Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue] = {}
    best_alignment_edge_type_map: Dict[str, str] = {}
    best_alignment_set_type_map: Dict[FrozenSet[CompositeEdge], str] = {}
    max_cost = 0

    alignment_sets: List[FrozenSet[CompositeEdge]] = []
    set_matches: Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue] = {}
    alignment_edge_type_map: Dict[str, str] = {}
    alignment_set_type_map: Dict[FrozenSet[CompositeEdge], str] = {}

    def is_feasible_extension(e1: CompositeEdge, e2: CompositeEdge, current_alignment: List[CompositeEdge],
                              sets_1: List[Set[Any]], sets_2: List[Set[Any]]) -> bool:
        # Check that the edge types match.
        if edge_type_map_1.get(str(e1)) != edge_type_map_2.get(str(e2)):
            logger.debug(f"Edge type mismatch: {e1} is {edge_type_map_1.get(str(e1))} vs {e2} is {edge_type_map_2.get(str(e2))}")
            return False

        if e1 is None or e2 is None:
            return False

        # Check compatibility
        incident_sets_1 = [s for s in sets_1 if e1 in s]
        incident_sets_2 = [s for s in sets_2 if e2 in s]

        if len(incident_sets_1) != 2 or len(incident_sets_2) != 2:
            logger.error("Error: Each edge must appear in exactly 2 sets (E-graph requirement)")
            return False

        # Check that the endpoint types match.
        types1 = frozenset(set_type_map_1.get(frozenset(s)) for s in incident_sets_1)
        types2 = frozenset(set_type_map_2.get(frozenset(s)) for s in incident_sets_2)
        logger.debug(f"types1={types1}")
        logger.debug(f"types2={types2}")
        if types1 != types2:
            logger.debug(f"Endpoint set type mismatch for edges {e1} and {e2}: {types1} vs {types2}")
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
            logger.debug(f"Neighbor mismatch: {aligned_neighbors_e1} vs {aligned_neighbors_e2}")
            return False

        x_y_in_set_matches = any(s == x for _, (s, t) in set_matches.items()) and \
                             any(s == y for _, (s, t) in set_matches.items())
        if x_y_in_set_matches:
            logger.debug("Edge already implied by existing set_matches for first graph.")
        a_b_in_set_matches = any(t == a for _, (s, t) in set_matches.items()) and \
                             any(t == b for _, (s, t) in set_matches.items())
        if a_b_in_set_matches:
            logger.debug("Edge already implied by existing set_matches for second graph.")

        if x_y_in_set_matches or a_b_in_set_matches:
            return False

        return True

    def backtrack_vf2(current_alignment: List[CompositeEdge],
                      frontier_edges_1: Set[CompositeEdge],
                      frontier_edges_2: Set[CompositeEdge],
                      matched_edges_1: Set[CompositeEdge],
                      matched_edges_2: Set[CompositeEdge]) -> None:
        nonlocal best_alignment, max_cost, best_set_matches, best_alignment_sets, best_alignment_edge_type_map, best_alignment_set_type_map
        nonlocal alignment_sets, set_matches, alignment_edge_type_map, alignment_set_type_map

        current_cost = sum(1 for edge in current_alignment if is_complete_edge(edge))
        if current_cost > max_cost:
            max_cost = current_cost
            best_alignment = list(current_alignment)
            best_alignment_sets = deepcopy(alignment_sets)
            best_set_matches = deepcopy(set_matches)
            best_alignment_edge_type_map = deepcopy(alignment_edge_type_map)
            best_alignment_set_type_map = deepcopy(alignment_set_type_map)

        #previous_edge_type_map = deepcopy(alignment_edge_type_map)

        for e1 in list(frontier_edges_1):
            for e2 in list(frontier_edges_2):
                logger.debug(f"EXPLORING candidate edge pair ({e1}, {e2})")
                if is_feasible_extension(e1, e2, current_alignment, sets_1, sets_2):
                    # Save current state before trying this candidate.
                    previous_current_alignment = deepcopy(current_alignment)
                    previous_alignment_sets = deepcopy(alignment_sets)
                    previous_set_matches = deepcopy(set_matches)
                    previous_edge_type_map = deepcopy(alignment_edge_type_map)
                    previous_set_type_map = deepcopy(alignment_set_type_map)

                    try:
                        logger.debug(f"Adding edge pair ({e1}, {e2}) to alignment...")
                        current_alignment.append((e1, e2))
                        alignment_edge_type_map[str((e1, e2))] = edge_type_map_1.get(str(e1))
                        # Update alignment sets and set matches.
                        alignment_sets, set_matches, alignment_edge_type_map, alignment_set_type_map = update_alignment_sets(
                            e1, e2, alignment_sets, set_matches, sets_1, sets_2,
                            current_alignment, edge_type_map_1, edge_type_map_2,
                            set_type_map_1, set_type_map_2, alignment_edge_type_map, alignment_set_type_map
                        )
                        # Check for implied edges. This is where the forbidden condition might be triggered.
                        set_matches, current_alignment, alignment_sets, alignment_edge_type_map, alignment_set_type_map = check_implied_edge_matches(
                            set_matches, current_alignment, alignment_sets,
                            edge_type_map_1, edge_type_map_2, set_type_map_1, set_type_map_2,
                            alignment_edge_type_map, alignment_set_type_map
                        )
                    except Exception as ex:
                        logger.debug(
                            f"Forbidden implied edge encountered when adding ({e1}, {e2}): {ex}. Rolling back and pruning this branch.")
                        # Restore previous state and skip further exploration of this branch.
                        current_alignment = previous_current_alignment
                        alignment_sets = previous_alignment_sets
                        set_matches = previous_set_matches
                        alignment_edge_type_map = previous_edge_type_map
                        alignment_set_type_map = previous_set_type_map
                        continue

                    alignment_set_type_map = {k: v for k, v in alignment_set_type_map.items() if k in alignment_sets}
                    logger.debug(f"Updated after adding ({e1}, {e2}):")
                    logger.debug(f"  Alignment Sets: {alignment_sets}")
                    logger.debug(f"  Set Matches: {set_matches}")
                    logger.debug(f"  Alignment Edge Types: {alignment_edge_type_map}")
                    logger.debug(f"  Alignment Set Types: {alignment_set_type_map}")
                    logger.debug(f"  Alignment edges: {current_alignment}")


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

                    logger.debug(f"Removing edge pair ({e1}, {e2}) from alignment, restoring state...")
                    # Restore state after backtracking from this candidate.
                    alignment_sets = previous_alignment_sets
                    set_matches = previous_set_matches
                    current_alignment = previous_current_alignment
                    alignment_edge_type_map = previous_edge_type_map
                    alignment_set_type_map = previous_set_type_map

                    matched_edges_1.remove(e1)
                    matched_edges_2.remove(e2)
                    frontier_edges_1.add(e1)
                    frontier_edges_2.add(e2)
        logger.debug(f"BEST alignment so far: {best_alignment}")
        logger.debug(f"best_alignment_edge_type_map: {best_alignment_edge_type_map}")

    frontier_edges_1 = set(edges_1)
    frontier_edges_2 = set(edges_2)
    logger.debug(f"Starting backtracking with initial frontier_edges_1: {frontier_edges_1}")
    backtrack_vf2([], frontier_edges_1, frontier_edges_2, set(), set())
    logger.debug(f"BEST alignment after backtracking: {best_alignment}")
    logger.debug(f"BEST set_matches after backtracking: {best_set_matches}")
    logger.debug(f"BEST alignment_edge_type_map: {best_alignment_edge_type_map}")
    logger.debug(f"BEST alignment_set_type_map: {best_alignment_set_type_map}")

    best_alignment, best_alignment_sets, best_set_matches, best_alignment_edge_type_map, best_alignment_set_type_map = handling_set_matches_isolated_edges(
        best_alignment, best_alignment_sets, best_set_matches, best_alignment_edge_type_map, best_alignment_set_type_map,
        sets_1, sets_2,
        edge_type_map_1, edge_type_map_2, set_type_map_1, set_type_map_2
    )

    best_alignment, best_alignment_sets, best_set_matches, best_alignment_edge_type_map, best_alignment_set_type_map = refine_alignment_with_unaligned_edges(
        best_alignment, best_alignment_sets, best_set_matches, best_alignment_edge_type_map, best_alignment_set_type_map, sets_1, sets_2,
        set(edges_1), set(edges_2), current_depth,
        edge_type_map_1, edge_type_map_2, set_type_map_1, set_type_map_2
    )

    logger.debug(f"Best alignment after refinement: {best_alignment}")
    logger.debug(f"Alignment Sets: {best_alignment_sets}")
    logger.debug(f"Cost: {max_cost}")
    logger.debug(f"Set Matches: {best_set_matches}")
    logger.debug(f"Alignment Edge Types: {best_alignment_edge_type_map}")
    logger.debug(f"Alignment Set Types: {best_alignment_set_type_map}")

    return best_alignment, best_alignment_sets, max_cost, best_set_matches, best_alignment_edge_type_map, best_alignment_set_type_map


def refine_alignment_with_unaligned_edges(
        alignment_edges: List[CompositeEdge],
        alignment_sets: List[FrozenSet[CompositeEdge]],
        set_matches: Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],
        alignment_edge_type_map: Dict[str, str],
        alignment_set_type_map: Dict[FrozenSet[CompositeEdge], str],
        sets_1: List[Set[Any]],
        sets_2: List[Set[Any]],
        edges_1: Set[CompositeEdge],
        edges_2: Set[CompositeEdge],
        current_depth: int,
        edge_type_map_1: Dict[CompositeEdge, str],
        edge_type_map_2: Dict[CompositeEdge, str],
        set_type_map_1: Dict[FrozenSet[Any], str],
        set_type_map_2: Dict[FrozenSet[Any], str]
) -> Tuple[
    List[CompositeEdge],
    List[FrozenSet[CompositeEdge]],
    Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],
    Dict[str, str],
    Dict[FrozenSet[CompositeEdge], str]
]:
    logger.debug("=== START PROCESSING UNALIGNED EDGES IN FIRST E-GRAPH ===")
    logger.debug(f"set_matches: {set_matches}")
    logger.debug(f"alignment_edges: {alignment_edges}")

    alignment_set_type_map = {k: v for k, v in alignment_set_type_map.items() if k in alignment_sets}

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

            type_ep1 = set_type_map_1.get(frozenset(endpoint_1))
            type_ep2 = set_type_map_1.get(frozenset(endpoint_2))

            w1_entry = next(
                ((s_key, t_key) for (s_key, t_key), (g, h) in set_matches.items() if g == frozenset(endpoint_1)), None)
            w2_entry = next(
                ((s_key, t_key) for (s_key, t_key), (g, h) in set_matches.items() if g == frozenset(endpoint_2)), None)
            # If found, unpack the alignment set part.
            w1 = w1_entry[0] if w1_entry else None
            w2 = w2_entry[0] if w2_entry else None
            logger.debug(f"w1 = {w1}\nw2 = {w2}")

            if w1 and w2:
                logger.debug(f"Adding ({e}, None) to two related sets.")
                updated_w1 = set(w1)
                updated_w2 = set(w2)
                updated_w1.add((e, None))
                updated_w2.add((e, None))
                new_w1 = frozenset(updated_w1)
                new_w2 = frozenset(updated_w2)
                alignment_sets[alignment_sets.index(w1)] = new_w1
                alignment_sets[alignment_sets.index(w2)] = new_w2
                alignment_set_type_map[new_w1] = type_ep1
                alignment_set_type_map[new_w2] = type_ep2

                if w1_entry is not None and (w1, w1_entry[1]) in set_matches:
                    old_value = set_matches.pop((w1, w1_entry[1]))
                    set_matches[(new_w1, type_ep1)] = old_value
                    logger.debug(
                        f"Updated set_matches: replaced key {(w1, w1_entry[1])} with {(new_w1, type_ep1)} -> {old_value}")
                else:
                    logger.error(f"No matching set found for endpoint_1: {endpoint_1}")

                if w2_entry is not None and (w2, w2_entry[1]) in set_matches:
                    old_value = set_matches.pop((w2, w2_entry[1]))
                    set_matches[(new_w2, type_ep2)] = old_value
                    logger.debug(
                        f"Updated set_matches: replaced key {(w2, w2_entry[1])} with {(new_w2, type_ep2)} -> {old_value}")
                else:
                    logger.error(f"No matching set found for endpoint_2: {endpoint_2}")

                alignment_edges.append((e, None))
                alignment_edge_type_map[str((e, None))] = edge_type_map_1[str(e)]
                U1.remove(e)
                progress = True

            elif w1 or w2:
                logger.debug(f"Adding ({e}, None) to one related set and creating a new set.")
                existing_entry = w1_entry if w1_entry else w2_entry
                existing_w = existing_entry[0]
                updated_w = set(existing_w)
                updated_w.add((e, None))
                new_w = frozenset(updated_w)
                alignment_sets[alignment_sets.index(existing_w)] = new_w
                alignment_set_type_map[new_w] = existing_entry[1]
                matching_entry = next(((s_key, t_key) for (s_key, t_key), (g, h) in set_matches.items()
                                       if g == frozenset(endpoint_1 if w1_entry else endpoint_2)), None)
                logger.debug(f"matching_entry = {matching_entry}")

                if matching_entry:
                    old_value = set_matches.pop(matching_entry)
                    set_matches[(new_w, matching_entry[1])] = old_value
                    logger.debug(
                        f"Updated set_matches: replaced key {matching_entry} with {(new_w, matching_entry[1])} -> {old_value}")
                else:
                    logger.error("No matching set found for endpoint_1 or endpoint_2.")

                new_other_w = frozenset({(e, None)})
                alignment_sets.append(new_other_w)

                if w1_entry:
                    alignment_set_type_map[new_other_w] = type_ep2
                    set_matches[(new_other_w, type_ep2)] = (frozenset(endpoint_2), 'nothing')
                    logger.debug(f"Created new set_matches entry: {new_other_w} -> ({endpoint_2}, 'nothing')")
                else:
                    alignment_set_type_map[new_other_w] = type_ep1
                    set_matches[(new_other_w, type_ep1)] = (frozenset(endpoint_1), 'nothing')
                    logger.debug(f"Created new set_matches entry: {new_other_w} -> ({endpoint_1}, 'nothing')")

                alignment_edges.append((e, None))
                alignment_edge_type_map[str((e, None))] = edge_type_map_1[str(e)]
                U1.remove(e)
                progress = True

            else:
                logger.debug(f"Could not find a related set for {e} in first E-graph. It is postponed.")
                continue

        if not progress:
            logger.error(f"Could not find a related set for any of {U1} in first E-graph.")
            break

    alignment_set_type_map = {k: v for k, v in alignment_set_type_map.items() if k in alignment_sets}

    logger.debug("=== START PROCESSING UNALIGNED EDGES IN SECOND E-GRAPH ===")
    logger.debug(f"alignment_sets = {alignment_sets}")
    logger.debug(f"alignment_edges = {alignment_edges}")
    logger.debug(f"set_matches = {set_matches}")

    aligned_edges_2 = {e2 for _, e2 in alignment_edges if e2 is not None}
    logger.debug(f"set(edges_2): {set(edges_2)}, aligned_edges_2: {aligned_edges_2}")
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

            type_ep1 = set_type_map_2.get(frozenset(endpoint_1))
            type_ep2 = set_type_map_2.get(frozenset(endpoint_2))

            w1_entry = next(
                ((s_key, t_key) for (s_key, t_key), (g, h) in set_matches.items() if h == frozenset(endpoint_1)), None)
            w2_entry = next(
                ((s_key, t_key) for (s_key, t_key), (g, h) in set_matches.items() if h == frozenset(endpoint_2)), None)

            # If found, unpack the alignment set part.
            w1 = w1_entry[0] if w1_entry else None
            w2 = w2_entry[0] if w2_entry else None
            logger.debug(f"w1 = {w1}\nw2 = {w2}")

            if w1 and w2:
                logger.debug(f"Adding (None, {e}) to two related sets.")
                updated_w1 = set(w1)
                updated_w2 = set(w2)
                new_edge = pad_edge(e, current_depth)
                updated_w1.add(new_edge)
                updated_w2.add(new_edge)
                new_w1 = frozenset(updated_w1)
                new_w2 = frozenset(updated_w2)
                alignment_sets[alignment_sets.index(w1)] = new_w1
                alignment_sets[alignment_sets.index(w2)] = new_w2
                alignment_set_type_map[new_w1] = type_ep1
                alignment_set_type_map[new_w2] = type_ep2

                if w1_entry is not None and (w1, w1_entry[1]) in set_matches:
                    old_value = set_matches.pop((w1, w1_entry[1]))
                    set_matches[(new_w1, type_ep1)] = old_value
                    logger.debug(
                        f"Updated set_matches: replaced key {(w1, w1_entry[1])} with {(new_w1, type_ep1)} -> {old_value}")
                else:
                    logger.error(f"No matching set found for endpoint_1: {endpoint_1}")

                if w2_entry is not None and (w2, w2_entry[1]) in set_matches:
                    old_value = set_matches.pop((w2, w2_entry[1]))
                    set_matches[(new_w2, type_ep2)] = old_value
                    logger.debug(
                        f"Updated set_matches: replaced key {(w2, w2_entry[1])} with {(new_w2, type_ep2)} -> {old_value}")
                else:
                    logger.error(f"No matching set found for endpoint_2: {endpoint_2}")

                new_edge = pad_edge(e, current_depth)
                alignment_edges.append(new_edge)
                alignment_edge_type_map[str(new_edge)] = edge_type_map_2[str(e)]
                U2.remove(e)
                progress = True
            elif w1 or w2:
                logger.debug(f"Adding ({e}, None) to one related set and creating a new set.")
                existing_entry = w1_entry if w1_entry else w2_entry
                existing_w = existing_entry[0]
                updated_w = set(existing_w)
                updated_w.add(pad_edge(e, current_depth))
                new_w = frozenset(updated_w)
                alignment_sets[alignment_sets.index(existing_w)] = new_w
                alignment_set_type_map[new_w] = existing_entry[1]
                matching_entry = next(((s_key, t_key) for (s_key, t_key), (g, h) in set_matches.items()
                                       if h == frozenset(endpoint_1 if w1_entry else endpoint_2)), None)
                logger.debug(f"matching_entry = {matching_entry}")
                if matching_entry:
                    old_value = set_matches.pop(matching_entry)
                    set_matches[(new_w, matching_entry[1])] = old_value
                    logger.debug(
                        f"Updated set_matches: replaced key {matching_entry} with {(new_w, matching_entry[1])} -> {old_value}")
                else:
                    logger.error("No matching set found for endpoint_1 or endpoint_2.")

                new_other_w = frozenset({pad_edge(e, current_depth)})
                alignment_sets.append(new_other_w)

                if w1_entry:
                    alignment_set_type_map[new_other_w] = type_ep2
                    set_matches[(new_other_w, type_ep2)] = ('nothing', frozenset(endpoint_2))
                    logger.debug(f"Created new set_matches entry: {new_other_w} -> ('nothing', {endpoint_2})")
                else:
                    alignment_set_type_map[new_other_w] = type_ep1
                    set_matches[(new_other_w, type_ep1)] = ('nothing', frozenset(endpoint_1))
                    logger.debug(f"Created new set_matches entry: {new_other_w} -> ('nothing', {endpoint_1})")

                new_edge = pad_edge(e, current_depth)
                alignment_edges.append(new_edge)
                alignment_edge_type_map[str(new_edge)] = edge_type_map_2[str(e)]
                U2.remove(e)
                progress = True

            else:
                logger.debug(f"Could not find a related set for {e} in second E-graph. It is postponed.")
                continue

        if not progress:
            logger.error(f"Could not find a related set for any of {U2} in second E-graph.")
            break

    alignment_set_type_map = {k: v for k, v in alignment_set_type_map.items() if k in alignment_sets}
    logger.debug("=== FINISHED PROCESSING UNALIGNED EDGES ===")
    logger.debug(f"Final alignment_edges: {alignment_edges}")
    logger.debug(f"Final alignment_sets: {alignment_sets}")
    logger.debug(f"Final set_matches: {set_matches}")

    return alignment_edges, alignment_sets, set_matches, alignment_edge_type_map, alignment_set_type_map


def align_multiple_e_graphs(
        edges_list: List[List[Dict[str, Any]]],
        sets_list: List[List[Dict[str, Any]]]
) -> Tuple[
    List[CompositeEdge],
    List[FrozenSet[CompositeEdge]],
    Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue],
    Dict[str, str],
    Dict[FrozenSet[CompositeEdge], str]
]:
    # Convert each e-graph's edge list (list of dicts) into a list of edge IDs and an edge type map.
    internal_edges_list = []
    edge_type_maps = []
    for e_graph in edges_list:
        edge_ids = []
        edge_type_map = {}
        for edge in e_graph:
            eid = edge["id"]
            edge_ids.append(eid)
            edge_type_map[eid] = edge.get("type")
        internal_edges_list.append(edge_ids)
        edge_type_maps.append(edge_type_map)

    # Convert each e-graph's sets (list of dicts) into a list of frozensets and a set type map.
    internal_sets_list = []
    set_type_maps = []
    for s_graph in sets_list:
        internal_sets = []
        set_type_map = {}
        for s in s_graph:
            edge_set = frozenset(s["edges"])
            internal_sets.append(edge_set)
            set_type_map[edge_set] = s.get("set_type")
        internal_sets_list.append(internal_sets)
        set_type_maps.append(set_type_map)

    # Check for duplicate edge labels across e-graphs.
    global_edge_to_graph: Dict[Any, int] = {}
    for i, graph_edges in enumerate(internal_edges_list):
        for edge in graph_edges:
            if edge in global_edge_to_graph:
                raise ValueError(
                    f"Edge label '{edge}' appears in both e-graph {global_edge_to_graph[edge]} and e-graph {i}. "
                    "It is recommended to use distinct edge labels."
                )
            global_edge_to_graph[edge] = i

    # For each e-graph, ensure each edge appears in exactly 0 or 2 sets.
    for i, (graph_edges, graph_sets) in enumerate(zip(internal_edges_list, internal_sets_list)):
        for edge in graph_edges:
            appearance_count = sum(1 for s in graph_sets if edge in s)
            if appearance_count == 0:
                logger.debug(f"In e-graph {i}, edge '{edge}' appears 0 times. Adding two singleton sets for it.")
                singleton = frozenset({edge})
                graph_sets.append(singleton)
                graph_sets.append(singleton)
                set_type_maps[i][singleton] = None
            elif appearance_count != 2:
                raise ValueError(
                    f"In e-graph {i}, edge '{edge}' appears {appearance_count} times. Each edge must appear in either 0 or 2 sets."
                )

    # Perform pairwise alignment.
    composite_edges = internal_edges_list[0]
    composite_sets = internal_sets_list[0]
    composite_edge_type_map = edge_type_maps[0]
    composite_set_type_map = set_type_maps[0]
    composite_set_matches: Dict[Tuple[FrozenSet[CompositeEdge], str], SetMatchValue] = {}
    current_depth = 1
    num_graphs = len(internal_edges_list)
    for i in range(1, num_graphs):
        current_depth += 1
        logger.info(f"Aligning e-graph number {current_depth}")
        logger.debug(f"composite_edge_type_map {composite_edge_type_map}")
        logger.debug(f"internal_edges_list[i]: {internal_edges_list[i]}")
        logger.debug(f"edge_type_maps[i]: {edge_type_maps[i]}")
        composite_edges, composite_sets, _, composite_set_matches, composite_edge_type_map, composite_set_type_map = find_best_alignment(
            composite_edges,
            composite_sets,
            internal_edges_list[i],
            internal_sets_list[i],
            current_depth,
            composite_edge_type_map,
            edge_type_maps[i],
            composite_set_type_map,
            set_type_maps[i]
        )
    return composite_edges, composite_sets, composite_set_matches, composite_edge_type_map, composite_set_type_map


def load_data_from_json(file_path: str) -> Tuple[List[List[Dict[str, Any]]], List[List[Dict[str, Any]]]]:
    with open(file_path, 'r') as f:
        data = json.load(f)

    edges_list = data.get("edges_list", [])
    sets_list = data.get("sets_list", [])

    if not isinstance(edges_list, list) or not all(isinstance(e_graph, list) for e_graph in edges_list):
        raise ValueError("edges_list should be a list of lists")
    if not isinstance(sets_list, list) or not all(isinstance(s_graph, list) for s_graph in sets_list):
        raise ValueError("sets_list should be a list of lists")

    for i, e_graph in enumerate(edges_list):
        for edge in e_graph:
            if not isinstance(edge, dict) or "id" not in edge or "type" not in edge:
                raise ValueError(
                    f"In e-graph {i}, each edge must be a dictionary with keys 'id' and 'type'. Found: {edge}"
                )

    for i, s_graph in enumerate(sets_list):
        for s in s_graph:
            if not isinstance(s, dict) or "edges" not in s or "set_type" not in s:
                raise ValueError(
                    f"In e-graph {i}, each set must be a dictionary with keys 'edges' and 'set_type'. Found: {s}"
                )
            if not isinstance(s["edges"], list):
                raise ValueError(
                    f"In e-graph {i}, the 'edges' value in a set must be a list. Found: {s['edges']}"
                )

    return edges_list, sets_list


def main() -> None:
    parser = argparse.ArgumentParser(description="Align multiple E-graphs with type constraints.")
    parser.add_argument("--log-level", type=str, default="info",
                        help="Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL). Default is INFO.")
    parser.add_argument("--log-file", type=str, default=None,
                        help="Optional file to output logs.")
    parser.add_argument("--input", type=str, default=None,
                        help="Path to JSON file containing 'edges_list' and 'sets_list'.")
    args = parser.parse_args()
    configure_logging(args.log_level, args.log_file)

    if args.input:
        try:
            edges_list, sets_list = load_data_from_json(args.input)
            logger.info(f"Loaded input from {args.input}")
        except Exception as e:
            logger.error(f"Failed to load input data: {e}")
            return
    else:

        edges_list = [
            [  # First E-graph
                {"id": "e1", "type": "alpha"},
                {"id": "e2", "type": "alpha"},
                {"id": "e3", "type": "beta"},
                {"id": "e4", "type": "beta"}
            ],
            [  # Second E-graph
                {"id": "f1", "type": "alpha"},
                {"id": "f2", "type": "alpha"},
                {"id": "f3", "type": "alpha"},
                {"id": "f4", "type": "alpha"},
                {"id": "f5", "type": "beta"},
                {"id": "f6", "type": "beta"},
                {"id": "f7", "type": "beta"}
            ],
            [  # Third E-graph
                {"id": "g1", "type": "alpha"},
                {"id": "g2", "type": "alpha"},
                {"id": "g3", "type": "beta"},
                {"id": "g4", "type": "beta"}
            ],

            [  # Fourth E-graph
                {"id": "k1", "type": "beta"},
                {"id": "k2", "type": "beta"},
                {"id": "k3", "type": "alpha"},
                {"id": "k4", "type": "alpha"}
            ]
        ]

        sets_list = [
            [  # Sets for first E-graph
                {"edges": ["e1", "e2"], "set_type": "setType1"},
                {"edges": ["e2", "e3"], "set_type": "setType1"},
                {"edges": ["e4", "e3"], "set_type": "setType1"},
                {"edges": ["e1", "e4"], "set_type": "setType2"}
            ],
            [  # Sets for second E-graph
                {"edges": ["f1", "f2", "f5"], "set_type": "setType1"},
                {"edges": ["f2", "f3"], "set_type": "setType1"},
                {"edges": ["f3", "f4"], "set_type": "setType1"},
                {"edges": ["f4", "f1"], "set_type": "setType1"},
                {"edges": ["f5", "f6"], "set_type": "setType1"},
                {"edges": ["f6", "f7"], "set_type": "setType2"},
                {"edges": ["f7"], "set_type": "setType2"}
            ],
            [  # Sets for third E-graph
                {"edges": ["g1"], "set_type": "setType2"},
                {"edges": ["g2", "g3"], "set_type": "setType1"},
                {"edges": ["g4", "g3"], "set_type": "setType1"},
                {"edges": ["g1", "g2"], "set_type": "setType1"},
                {"edges": ["g4"], "set_type": "setType2"}
            ],


            [  # Sets for fourth E-graph
                {"edges": ["k1"], "set_type": "setType1"},
                {"edges": ["k2", "k3"], "set_type": "setType1"},
                {"edges": ["k4", "k3"], "set_type": "setType1"},
                {"edges": ["k1", "k2"], "set_type": "setType2"},
                {"edges": ["k4"], "set_type": "setType2"}
            ]
        ]

        logger.info("No input file provided. Using default test data.")

    try:
        alignment_edges, alignment_sets, set_matches, edge_types, set_types = align_multiple_e_graphs(edges_list, sets_list)
    except Exception as e:
        logger.error(f"Alignment failed: {e}")
        return

    logger.info("Final Alignment Edges:")
    for edge in alignment_edges:
        logger.info(edge)

    logger.info("Final Alignment Sets:")
    for idx, edge_set in enumerate(alignment_sets):
        logger.info(f"Set {idx + 1}: {edge_set}")

    logger.info("Final Set Matches:")
    for s3, (s1, s2) in set_matches.items():
        logger.info(f"{s3} -> {s1}, {s2}")

    logger.info("Final Edge types:")
    for key, value in edge_types.items():
        logger.info(f"  {key}: {value}")

    logger.info("Final Set types:")
    for key, value in set_types.items():
        logger.info(f"  {set(key)}: {value}")


if __name__ == "__main__":
    main()
