from collections import defaultdict, deque, Counter
from itertools import combinations
from copy import deepcopy
import networkx as nx
import matplotlib.pyplot as plt
import itertools
import argparse



def pad_none(depth):
    """
    Recursively produces a nested tuple of Nones of the given depth.
    For example, pad_none(1) returns None,
                 pad_none(2) returns (None, None),
                 pad_none(3) returns ((None, None), None), etc.
    """
    if depth <= 1:
        return None
    else:
        return (pad_none(depth - 1), None)


def pad_edge(e, current_depth):
    """
    Given a new edge 'e' and the current depth (i.e. the number of input e-graphs
    you are aligning), returns a nested tuple that pads the left side with the appropriate
    number of Nones.

    For example, if current_depth == 2, returns (None, e);
                 if current_depth == 3, returns ((None, None), e);
                 if current_depth == 4, returns (((None, None), None), e).
    """
    if current_depth <= 1:
        return e
    else:
        return (pad_none(current_depth - 1), e)


def refine_alignment_with_unaligned_edges(alignment_edges, alignment_sets, set_matches, sets_1, sets_2, edges_1, edges_2, current_depth):
    """
    Ensures that unaligned edges from both E-graphs are correctly added to alignment_edges and alignment_sets using set-matches.

    :param alignment_edges: List of already aligned edges (contains (e1, e2) pairs)
    :param alignment_sets: List of sets currently in the alignment
    :param set_matches: Dictionary tracking set matches
    :param sets_1: List of sets defining the first E-graph
    :param sets_2: List of sets defining the second E-graph
    :param edges_1: Set of edges in the first E-graph
    :param edges_2: Set of edges in the second E-graph
    :return: Updated alignment_edges, alignment_sets, set_matches
    """

    print("\n[DEBUG] --- START PROCESSING UNALIGNED EDGES IN FIRST E-GRAPH ---")
    print(f"set_matches+={set_matches}")
    print(f"alignment_edges+={alignment_edges}")

    # Step 1: Extract the set of unaligned edges in the first graph
    aligned_edges_1 = {e1 for e1, e2 in alignment_edges}
    U1 = set(edges_1) - aligned_edges_1  # Unaligned edges in the first graph

    print(f"[DEBUG] Initial unaligned edges in first E-graph (U1): {U1}")

    # Step 2: Process unaligned edges in the first E-graph
    while U1:
        print(f"U1={U1}")
        progress = False  # Flag to check if at least one valid edge is processed

        for e in list(U1):  # Iterate over a list copy to allow removal
            endpoint_sets = [s for s in sets_1 if e in s]
            print(f"endpoint_sets = {endpoint_sets}")
            print("Current set_matches:", set_matches)

            if len(endpoint_sets) != 2:
                print(f"[ERROR] Edge {e} should appear in exactly two sets, but found {len(endpoint_sets)}")
                continue

            endpoint_1, endpoint_2 = endpoint_sets

            # For first E-graph, match by comparing the first coordinate (the first element of the tuple)
            w1 = next((key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_1)), None)
            w2 = next((key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_2)), None)

            print(f"w1 = {w1}\n w2 = {w2}")


            if w1 and w2:
                print(f"[DEBUG] Adding ({e}, None) to two related sets.")

                updated_w1, updated_w2 = set(w1), set(w2)
                updated_w1.add((e, None))
                updated_w2.add((e, None))

                new_w1, new_w2 = frozenset(updated_w1), frozenset(updated_w2)
                alignment_sets[alignment_sets.index(w1)] = new_w1
                alignment_sets[alignment_sets.index(w2)] = new_w2

                # Update set_matches for both endpoints
                matching_entry_1 = next((key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_1)), None)
                matching_entry_2 = next((key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_2)), None)

                if matching_entry_1:
                    old_value = set_matches[matching_entry_1]
                    del set_matches[matching_entry_1]
                    set_matches[new_w1] = old_value
                    print(f"[DEBUG] Updated set_matches: Removed {matching_entry_1}, Added {new_w1} -> {old_value}")
                else:
                    print(f"[ERROR] No matching set found for endpoint_1: {endpoint_1}")

                if matching_entry_2:
                    old_value = set_matches[matching_entry_2]
                    del set_matches[matching_entry_2]
                    set_matches[new_w2] = old_value
                    print(f"[DEBUG] Updated set_matches: Removed {matching_entry_2}, Added {new_w2} -> {old_value}")
                else:
                    print(f"[ERROR] No matching set found for endpoint_2: {endpoint_2}")

                alignment_edges.append((e, None))
                U1.remove(e)
                progress = True

            elif w1 or w2:
                print(f"[DEBUG] Adding ({e}, None) to one related set and creating a new set.")

                existing_w = w1 if w1 else w2
                updated_w = set(existing_w)
                updated_w.add((e, None))
                new_w = frozenset(updated_w)
                alignment_sets[alignment_sets.index(existing_w)] = new_w

                # Update the matching entry for the aligned endpoint
                matching_entry = next(
                    (key for key, (s, w) in set_matches.items() if s == frozenset(endpoint_1 if w1 else endpoint_2)),
                    None
                )
                print(f"matching_entry = {matching_entry}")
                if matching_entry:
                    old_value = set_matches[matching_entry]
                    del set_matches[matching_entry]
                    set_matches[new_w] = old_value
                    print(f"[DEBUG] Updated set_matches: Removed {matching_entry}, Added {new_w} -> {old_value}")
                else:
                    print(f"[ERROR] No matching set found for endpoint_1 or endpoint_2.")

                # Create a new set for the missing endpoint
                new_other_w = frozenset({(e, None)})
                alignment_sets.append(new_other_w)


                # Instead of searching for a matching entry, create one directly:
                if w1:
                    # w1 exists: endpoint_1 is already aligned; so create a new entry using endpoint_1 as first coordinate.
                    set_matches[new_other_w] = (endpoint_2, "nothing")
                    print(f"[DEBUG] Created new set_matches entry: {new_other_w} -> ({endpoint_2}, 'nothing')")
                else:
                    # w2 exists: endpoint_2 is already aligned; so create a new entry using endpoint_2 as second coordinate.
                    set_matches[new_other_w] = (endpoint_1, "nothing")
                    print(f"[DEBUG] Created new set_matches entry: {new_other_w} -> ({endpoint_1}, 'nothing')")

                alignment_edges.append((e, None))
                U1.remove(e)
                progress = True

            else:
                print(f"[DEBUG] Could not find a related set for {e} in first E-graph. It is it is postponed.")
                # Do not mark progress, so that the edge remains in U1 for a later pass.
                continue

        # End of round for U1. If no edge was processed in this round, then we have a deadlock.
        if not progress:
            print(f"[ERROR] Could not find a related set for any of {U1} in first E-graph.")
            break

    print("\n[DEBUG] --- START PROCESSING UNALIGNED EDGES IN SECOND E-GRAPH ---")
    print(f"alignment_sets = {alignment_sets}")
    print(f"alignment_edges = {alignment_edges}")
    print(f"set_matches = {set_matches}")

    # Step 3: Extract the set of unaligned edges in the second graph
    aligned_edges_2 = {e2 for e1, e2 in alignment_edges}
    U2 = set(edges_2) - aligned_edges_2  # Unaligned edges in the second graph

    print(f"[DEBUG] Initial unaligned edges in second E-graph (U2): {U2}")

    # Step 4: Process unaligned edges in the second E-graph
    while U2:
        print(f"U2={U2}")
        progress = False

        for e in list(U2):
            endpoint_sets = [s for s in sets_2 if e in s]
            print(f"endpoint_sets = {endpoint_sets}")
            print("Current set_matches:", set_matches)

            if len(endpoint_sets) != 2:
                print(f"[ERROR] Edge {e} should appear in exactly two sets, but found {len(endpoint_sets)}")
                continue

            endpoint_1, endpoint_2 = endpoint_sets

            # For the second E-graph, match by comparing the second coordinate
            w1 = next((key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_1)), None)
            w2 = next((key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_2)), None)
            print(f"w1 = {w1}\n w2 = {w2}")

            if w1 and w2:
                print(f"[DEBUG] Adding (None, {e}) to two related sets.")

                updated_w1, updated_w2 = set(w1), set(w2)
                updated_w1.add(pad_edge(e, current_depth))
                updated_w2.add(pad_edge(e, current_depth))

                new_w1, new_w2 = frozenset(updated_w1), frozenset(updated_w2)
                alignment_sets[alignment_sets.index(w1)] = new_w1
                alignment_sets[alignment_sets.index(w2)] = new_w2

                # Update set_matches using the second coordinate (w) for matching
                matching_entry_1 = next((key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_1)),
                                        None)
                matching_entry_2 = next((key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_2)),
                                        None)
                if matching_entry_1:
                    old_value = set_matches[matching_entry_1]
                    del set_matches[matching_entry_1]
                    set_matches[new_w1] = old_value
                    print(f"[DEBUG] Updated set_matches: Removed {matching_entry_1}, Added {new_w1} -> {old_value}")
                else:
                    print(f"[ERROR] No matching set found for endpoint_1: {endpoint_1}")

                if matching_entry_2:
                    old_value = set_matches[matching_entry_2]
                    del set_matches[matching_entry_2]
                    set_matches[new_w2] = old_value
                    print(f"[DEBUG] Updated set_matches: Removed {matching_entry_2}, Added {new_w2} -> {old_value}")
                else:
                    print(f"[ERROR] No matching set found for endpoint_2: {endpoint_2}")

                alignment_edges.append(pad_edge(e, current_depth))
                U2.remove(e)
                progress = True

            elif w1 or w2:
                print(f"[DEBUG] Adding (None, {e}) to one related set and creating a new set.")

                existing_w = w1 if w1 else w2
                updated_w = set(existing_w)
                updated_w.add(pad_edge(e, current_depth))
                new_w = frozenset(updated_w)
                alignment_sets[alignment_sets.index(existing_w)] = new_w

                # Update the matching entry for the aligned endpoint
                matching_entry = next(
                    (key for key, (s, w) in set_matches.items() if w == frozenset(endpoint_1 if w1 else endpoint_2)),
                    None
                )
                print(f"matching_entry = {matching_entry}")
                if matching_entry:
                    old_value = set_matches[matching_entry]
                    del set_matches[matching_entry]
                    set_matches[new_w] = old_value
                    print(f"[DEBUG] Updated set_matches: Removed {matching_entry}, Added {new_w} -> {old_value}")
                else:
                    print(f"[ERROR] No matching set found for endpoint_2 or endpoint_1.")

                # Create a new set for the missing endpoint
                new_other_w = frozenset({pad_edge(e, current_depth)})
                alignment_sets.append(new_other_w)

                # For the second E-graph, if one endpoint is matched then create a new entry using the missing coordinate.
                if w1:
                    # w1 exists: endpoint_1 is already aligned; so the missing second coordinate comes from endpoint_2.
                    set_matches[new_other_w] = ("nothing", endpoint_2)
                    print(f"[DEBUG] Created new set_matches entry: {new_other_w} -> ('nothing', {endpoint_2})")
                else:
                    # w2 exists: endpoint_2 is already aligned; so the missing second coordinate comes from endpoint_1.
                    set_matches[new_other_w] = ("nothing", endpoint_1)
                    print(f"[DEBUG] Created new set_matches entry: {new_other_w} -> ('nothing', {endpoint_1})")

                alignment_edges.append(pad_edge(e, current_depth))
                U2.remove(e)
                progress = True

            else:
                print(f"[DEBUG] Could not find a related set for {e} in second E-graph. It is it is postponed.")
                continue

        if not progress:
            print(f"[ERROR] Could not find a related set for any of {U2} in second E-graph.")
            break

    print("\n[DEBUG] --- FINISHED PROCESSING UNALIGNED EDGES ---")
    print(f"[DEBUG] Final alignment_edges: {alignment_edges}")
    print(f"[DEBUG] Final alignment_sets: {alignment_sets}")
    print(f"[DEBUG] Final set_matches: {set_matches}")

    return alignment_edges, alignment_sets, set_matches


# def check_implied_edge_matches(set_matches, alignment_edges, alignment_sets):
#     print("\n[DEBUG] Checking for implied edge matches...")
#
#     for (w1, (v1, u1)), (w2, (v2, u2)) in combinations(set_matches.items(), 2):
#         # Ensure we are comparing two distinct alignment sets
#         if w1 != w2:
#             # Check for intersection between v1, v2 and u1, u2
#             v1_set, u1_set = frozenset(v1), frozenset(u1)
#             v2_set, u2_set = frozenset(v2), frozenset(u2)
#
#             if v1_set.intersection(v2_set) and u1_set.intersection(u2_set):
#                 # Extract the implied edge
#                 e1 = list(v1_set.intersection(v2_set))[0]
#                 e2 = list(u1_set.intersection(u2_set))[0]
#
#                 # Check if the implied edge is already in alignment_sets
#                 if not any((e1, e2) in s for s in alignment_sets):
#                     print(f"  [DEBUG] Found implied edge: ({e1}, {e2})")
#                     print(set_matches)
#                     alignment_edges.append((e1, e2))
#
#                     # Add the implied edge to both w1 and w2
#                     updated_w1 = set(w1)
#                     updated_w1.add((e1, e2))
#                     alignment_sets[alignment_sets.index(w1)] = frozenset(updated_w1)
#                     #set_matches[frozenset(updated_w1)] = (v1, u1)
#
#                     updated_w2 = set(w2)
#                     updated_w2.add((e1, e2))
#                     alignment_sets[alignment_sets.index(w2)] = frozenset(updated_w2)
#                     #set_matches[frozenset(updated_w2)] = (v2, u2)
#
#                     # Explicitly update set_matches for w1
#                     if frozenset(w1) in set_matches:
#                         set_matches[frozenset(updated_w1)] = set_matches.pop(
#                             frozenset(w1))  # Replace key with updated set
#
#                     # Explicitly update set_matches for w2
#                     if frozenset(w2) in set_matches:
#                         set_matches[frozenset(updated_w2)] = set_matches.pop(
#                             frozenset(w2))  # Replace key with updated set
#
#                     print(f"  [DEBUG] Updated w1: {frozenset(updated_w1)}")
#                     print(f"  [DEBUG] Updated w2: {frozenset(updated_w2)}")
#
#     print(f"[DEBUG] Updated alignment edges: {alignment_edges}")
#     print(f"[DEBUG] Updated alignment sets: {alignment_sets}")
#     print(f"[DEBUG] Updated set matches: {set_matches}")
#     return set_matches, alignment_edges, alignment_sets

def check_implied_edge_matches(set_matches, alignment_edges, alignment_sets):
    print("\n[DEBUG] Checking for implied edge matches...")

    # Continue looping until no updates occur in a full pass.
    updated = True
    while updated:
        updated = False  # Assume no updates, then check each pair.

        # Iterate over all distinct pairs in set_matches.
        for (w1, (v1, u1)), (w2, (v2, u2)) in combinations(set_matches.items(), 2):
            if w1 != w2:
                # Convert endpoint sets to frozensets (if they aren't already)
                v1_set, u1_set = frozenset(v1), frozenset(u1)
                v2_set, u2_set = frozenset(v2), frozenset(u2)

                # If both intersections are non-empty, then an implied edge exists.
                if v1_set.intersection(v2_set) and u1_set.intersection(u2_set):
                    e1 = list(v1_set.intersection(v2_set))[0]
                    e2 = list(u1_set.intersection(u2_set))[0]

                    # If this implied edge isn't already present, then update.
                    if not any((e1, e2) in s for s in alignment_sets):
                        print(f"  [DEBUG] Found implied edge: ({e1}, {e2})")
                        print("         set_matches =", set_matches)
                        alignment_edges.append((e1, e2))

                        # Update w1 and w2 by taking the union with the new edge.
                        updated_w1 = set(w1)
                        updated_w1.add((e1, e2))
                        new_w1 = frozenset(updated_w1)

                        updated_w2 = set(w2)
                        updated_w2.add((e1, e2))
                        new_w2 = frozenset(updated_w2)

                        # Update alignment_sets: replace the old sets with the new ones.
                        alignment_sets[alignment_sets.index(w1)] = new_w1
                        alignment_sets[alignment_sets.index(w2)] = new_w2

                        # Update set_matches keys by re-keying.
                        if frozenset(w1) in set_matches:
                            set_matches[frozenset(new_w1)] = set_matches.pop(frozenset(w1))
                        if frozenset(w2) in set_matches:
                            set_matches[frozenset(new_w2)] = set_matches.pop(frozenset(w2))

                        print(f"  [DEBUG] Updated w1: {new_w1}")
                        print(f"  [DEBUG] Updated w2: {new_w2}")

                        # Mark that we did an update.
                        updated = True
                        break  # Break out of the for-loop immediately to avoid the key finding issue (that key might no longer find it.)

        if updated:
            print("[DEBUG] One round complete. Restarting iteration due to updates.")
        else:
            print("[DEBUG] No new implied edge found in this round. Terminating.")

    print(f"[DEBUG] Final alignment edges: {alignment_edges}")
    print(f"[DEBUG] Final alignment sets: {alignment_sets}")
    print(f"[DEBUG] Final set matches: {set_matches}")
    return set_matches, alignment_edges, alignment_sets


def update_alignment_sets(e1, e2, alignment_sets, set_matches, sets_1, sets_2, isolated_edges, current_alignment):
    print(f"\n[DEBUG] Updating alignment sets for edge ({e1}, {e2})...")

    # Step 1: Check for Isolation and Identify (f1, f2)
    f1, f2 = None, None
    is_isolated = True

    for alignment_set in alignment_sets:
        for pair in alignment_set:
            # Explicitly check S1 and S2 for connectivity
            if any(e1 in s and pair[0] in s for s in sets_1) and any(e2 in s and pair[1] in s for s in sets_2):
                f1, f2 = pair
                is_isolated = False
                break
        if not is_isolated:
            break

    if is_isolated:
        print("  [DEBUG] Edge is isolated. Deferring decision.")
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
        print(f"  [DEBUG] Isolated edge processed. Updated alignment sets: {alignment_sets}")
        return alignment_sets, set_matches  # Exit early for isolated edges

    # Ensure that we have a valid (f1, f2) for a non-isolated edge
    assert f1 is not None and f2 is not None, "[ERROR] Non-isolated edge but no (f1, f2) found!"
    print(f"[DEBUG] Fixed aligned edge (f1, f2): ({f1}, {f2})")

    # Step 2: Define v1, v2, v3, u1, u2, u3
    v1 = frozenset(next((s for s in sets_1 if e1 in s and f1 not in s), set()))
    v2 = frozenset(next((s for s in sets_1 if e1 in s and f1 in s), set()))
    v3 = frozenset(next((s for s in sets_1 if f1 in s and e1 not in s), set()))

    u1 = frozenset(next((s for s in sets_2 if e2 in s and f2 not in s), set()))
    u2 = frozenset(next((s for s in sets_2 if e2 in s and f2 in s), set()))
    u3 = frozenset(next((s for s in sets_2 if f2 in s and e2 not in s), set()))

    print(f"  [DEBUG] Unique sets for (f1, f2): v1={v1}, v2={v2}, v3={v3}, u1={u1}, u2={u2}, u3={u3}")

    # Step 3: Process Alignment Updates
    set_match_v1 = any(v1 in match for match in set_matches.values())
    set_match_v2 = any(v2 in match for match in set_matches.values())
    set_match_v3 = any(v3 in match for match in set_matches.values())

    print(f"  [DEBUG] Current set matches: {set_matches}")
    print(f"  [DEBUG] Current alignment sets: {alignment_sets}")

    if not set_match_v1 and not set_match_v2 and not set_match_v3:
        print("  [DEBUG] Case A: No set-match for v1, v2, or v3.")
        # Handle Case A
        # Find all sets in alignment_sets that are equal to {(f1, f2)}
        singleton_sets = [s for s in alignment_sets if s == frozenset({(f1, f2)})]

        # Ensure there are at least two such sets
        if len(singleton_sets) >= 2:
            w2 = singleton_sets[0]
            w3 = singleton_sets[1]

            print(f"  [DEBUG] Found two equal sets for (f1, f2): w2={w2}, w3={w3}")

            # Add (e1, e2) to w2 only (not w3)
            updated_w2 = set(w2)
            updated_w2.add((e1, e2))
            alignment_sets[alignment_sets.index(w2)] = frozenset(updated_w2)

            print(f"  [DEBUG] Updated w2: {frozenset(updated_w2)} (added (e1, e2))")
            print(f"  [DEBUG] w3 remains unchanged: {w3}")

            # Update set-matches for w2 and w3
            set_matches[frozenset(updated_w2)] = (v2, u2)  # v2 to u2 corresponds to w2
            set_matches[w3] = (v3, u3)  # v3 to u3 corresponds to w3

            print(f"  [DEBUG] Updated set-matches: w2 -> (v2, u2), w3 -> (v3, u3)")
        else:
            print(f"[ERROR] Could not find two equal sets for (f1, f2) in alignment_sets.")

        # Step 5: If v1 intersects with already matched edges
        v1_aligned_edges_count = len({e for e in v1 if any(e == edge[0] for edge in current_alignment)})
        print(f"v1_aligned_edges_count={v1_aligned_edges_count}")
        if v1_aligned_edges_count > 2:
            print("Error: There has to be a matching as there are two edges on e1 in the alignment")
        elif v1_aligned_edges_count == 1:
            print(f"current_alignment={current_alignment}")
            w1 = frozenset({(e1, e2)})
            alignment_sets.append(w1)
            set_matches[w1] = (v1, u1)
            print(f"  [DEBUG] Created new set w1: {w1}")
        elif v1_aligned_edges_count == 2:
            print("g1g2 Needed!!!!!")
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
                    print(f"  [DEBUG] Updated w1: {frozenset(updated_w1)} (added (e1, e2))")
                    print(f"  [DEBUG] w4 remains unchanged: {w4}")
                    v0 = frozenset(next((s for s in sets_1 if g1 in s and e1 not in s), set()))
                    u0 = frozenset(next((s for s in sets_2 if g2 in s and e2 not in s), set()))
                    set_matches[frozenset(updated_w1)] = (v1, u1)
                    set_matches[w4] = (v0, u0)
                    print(f"  [DEBUG] Updated set-matches: w1 -> (v1, u1), w4 -> (v0, u0)")

        print(f"[DEBUG] Final alignmentAA: {current_alignment}")
        print(f"[DEBUG] Final alignment sets: {alignment_sets}")
        print(f"[DEBUG] Final set matches: {set_matches}")
        return alignment_sets, set_matches
    else:
        print("  [DEBUG] Case B: Found set-match for v1, v2, or v3.")
        # Handle subcases within Case B

        # Subcase B.1: Set-match includes v1
        if set_match_v1:
            w1 = next((s3 for s3, (s1, _) in set_matches.items() if s1 == v1), None)
            if w1:
                updated_w1 = set(w1)
                updated_w1.add((e1, e2))
                alignment_sets[alignment_sets.index(w1)] = frozenset(updated_w1)

                # Ensure set_matches is updated for w1
                set_matches[frozenset(updated_w1)] = set_matches.pop(w1)

                print(f"  [DEBUG] Updated w1: {frozenset(updated_w1)} (added (e1, e2))")

            # Update w2 and w3
            singleton_sets = [s for s in alignment_sets if s == frozenset({(f1, f2)})]
            if len(singleton_sets) >= 2:
                w2 = singleton_sets[0]
                w3 = singleton_sets[1]
                updated_w2 = set(w2)
                updated_w2.add((e1, e2))
                alignment_sets[alignment_sets.index(w2)] = frozenset(updated_w2)
                set_matches[frozenset(updated_w2)] = (v2, u2)
                set_matches[w3] = (v3, u3)
                print(f"  [DEBUG] Updated w2: {frozenset(updated_w2)}")
                print(f"  [DEBUG] w3 remains unchanged: {w3}")

        # Subcase B.2: Set-match includes v2
        elif set_match_v2:
            v1_aligned_edges = {e for e in v1 if any(e == edge[0] for edge in current_alignment)}
            print(f"v1_aligned_edges(Case B)={v1_aligned_edges}")
            print(f"current_alignment={current_alignment}")

            w2 = next((s3 for s3, (s1, _) in set_matches.items() if s1 == v2), None)
            if w2:
                updated_w2 = set(w2)
                updated_w2.add((e1, e2))
                alignment_sets[alignment_sets.index(w2)] = frozenset(updated_w2)
                # Ensure set_matches is updated for w2, as well
                set_matches[frozenset(updated_w2)] = set_matches.pop(w2)
                print(f"  [DEBUG] Updated w2: {frozenset(updated_w2)} (added (e1, e2))")

            # Handle |v1 ∩ already matched edges| conditions

            if len(v1_aligned_edges) == 1:
                w1 = frozenset({(e1, e2)})
                alignment_sets.append(w1)
                set_matches[w1] = (v1, u1)
                print(f"  [DEBUG] Created new set w1: {w1}")
            elif len(v1_aligned_edges) == 2:
                # Handle g1, g2 logic for w1 and w4

                g1 = next((e for e in v1 if e != e1 and any(e == edge[0] for edge in current_alignment)), None)
                g2 = next((e for e in u1 if e != e2 and any(e == edge[1] for edge in current_alignment)), None)

                print(f" We have g1={g1},g2={g2}")
                if g1 and g2:
                    equal_sets = [s for s in alignment_sets if s == frozenset({(g1, g2)})]
                    print(f"equal_sets={equal_sets}")
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
                        print(f"  [DEBUG] Updated w1: {frozenset(updated_w1)}")
                        print(f"  [DEBUG] Updated w4: {w4}")

    print(f"[DEBUG] Final alignmentAA: {current_alignment}")
    print(f"[DEBUG] Final alignment sets: {alignment_sets}")
    print(f"[DEBUG] Final set matches: {set_matches}")
    return alignment_sets, set_matches

def is_complete_edge(edge):
    """
    Recursively checks if an edge (which may be nested) contains any None values.
    Returns True if no None is found; otherwise, returns False.
    """
    # If the edge is a tuple or list, check each sub-element recursively.
    if isinstance(edge, (tuple, list)):
        return all(is_complete_edge(sub) for sub in edge)
    # Otherwise, just check if the element is not None.
    return edge is not None


def find_best_alignment(edges_1, sets_1, edges_2, sets_2, current_depth):
    best_alignment = []
    max_cost = 0
    best_alignment_sets = []
    best_set_matches = {}
    alignment_sets = []
    set_matches = {}
    isolated_edges = {}

    def is_feasible_extension(e1, e2, current_alignment, sets_1, sets_2):
        if e1 is None or e2 is None:
            return False

        # Step 1: Identify the two sets containing e1 in S1 and e2 in S2
        incident_sets_1 = [s for s in sets_1 if e1 in s]
        incident_sets_2 = [s for s in sets_2 if e2 in s]

        # Ensure each edge appears in exactly two sets
        if len(incident_sets_1) != 2 or len(incident_sets_2) != 2:
            print("Error: Inputs have to be E-graphs")
            return False

        x, y = incident_sets_1  # The two sets containing e1 in S1
        a, b = incident_sets_2  # The two sets containing e2 in S2

        # Step 2: Compute the union of neighbors
        union_neighbors_e1 = set(x).union(set(y)) - {e1}  # All edges in x ∪ y except e1
        union_neighbors_e2 = set(a).union(set(b)) - {e2}  # All edges in a ∪ b except e2

        # Step 3: Identify aligned neighbors
        aligned_neighbors_e1 = {
            (aligned_e1, aligned_e2) for aligned_e1, aligned_e2 in current_alignment if aligned_e1 in union_neighbors_e1
        }
        aligned_neighbors_e2 = {
            (aligned_e1, aligned_e2) for aligned_e1, aligned_e2 in current_alignment if aligned_e2 in union_neighbors_e2
        }
        print(f"aligned_neighbors_e1:{aligned_neighbors_e1} \n aligned_neighbors_e2:{aligned_neighbors_e2}")
        # Step 4: Compare aligned neighbors
        if aligned_neighbors_e1 != aligned_neighbors_e2:
            print(f"aligned_neighbors_e1:{aligned_neighbors_e1} is NOT = aligned_neighbors_e2:{aligned_neighbors_e2}")
            return False  # Not feasible if the aligned neighbors are inconsistent

        # Step 5: Check if both x and y are already in set_matches
        x_y_in_set_matches = any(s == x for _, (s, t) in set_matches.items()) and \
                             any(s == y for _, (s, t) in set_matches.items())

        if x_y_in_set_matches:
            print("Should be already there as an implied edge (for the first coordinate)")
        # Step 6: Check if both a and b are already in set_matches
        a_b_in_set_matches = any(t == a for _, (s, t) in set_matches.items()) and \
                             any(t == b for _, (s, t) in set_matches.items())

        if a_b_in_set_matches:
            print("Should be already there as an implied edge (for the second coordinate)")
        # If both x and y or both a and b are already in set_matches, return False
        if x_y_in_set_matches or a_b_in_set_matches:
            return False

        return True


    def backtrack_vf2(current_alignment, frontier_edges_1, frontier_edges_2, matched_edges_1, matched_edges_2):
        nonlocal best_alignment, max_cost, best_set_matches, best_alignment_sets, alignment_sets, set_matches, isolated_edges

        current_cost = sum(1 for edge in current_alignment if is_complete_edge(edge))


        if current_cost > max_cost:
            max_cost = current_cost
            best_alignment = list(current_alignment)
            best_alignment_sets = deepcopy(alignment_sets)
            best_set_matches = deepcopy(set_matches)

        for e1 in list(frontier_edges_1):
            for e2 in list(frontier_edges_2):
                print(f"\n[DEBUG] EXPLORING candidate edge pair ({e1}, {e2})")
                if is_feasible_extension(e1, e2, current_alignment, sets_1, sets_2):

                    # Save the current state
                    previous_current_alignment = deepcopy(current_alignment)

                    print(f"\n[DEBUG] Adding edge pair ({e1}, {e2}) to alignment...")
                    current_alignment.append((e1, e2))

                    # Save the current state
                    previous_alignment_sets = deepcopy(alignment_sets)
                    previous_set_matches = deepcopy(set_matches)


                    print(f"  [DEBUG] Saved state before adding ({e1}, {e2}):")
                    print(f"    Alignment Sets: {previous_alignment_sets}")
                    print(f"    Set Matches: {previous_set_matches}")
                    print(f"    previous_current_alignment: {previous_current_alignment}")

                    # Update alignment sets and set matches
                    alignment_sets, set_matches = update_alignment_sets(
                        e1, e2, alignment_sets, set_matches, sets_1, sets_2, isolated_edges, current_alignment
                    )

                    set_matches, current_alignment, alignment_sets = check_implied_edge_matches(set_matches,
                                                                                              current_alignment,
                                                                                              alignment_sets)

                    print(f"  [DEBUG] Updated after adding ({e1}, {e2}):")
                    print(f"    Alignment Sets: {alignment_sets}")
                    print(f"    Set Matches: {set_matches}")

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

                    # Restore the state
                    print(f"\n[DEBUG] Removing edge pair ({e1}, {e2}) from alignment...")
                    alignment_sets = previous_alignment_sets
                    set_matches = previous_set_matches
                    print(f"current_alignment_before POP={current_alignment}")
                    #current_alignment.pop()
                    current_alignment = previous_current_alignment
                    print(f"current_alignment_after POP={current_alignment}")


                    print(f"  [DEBUG] Restored state after removing ({e1}, {e2}):")
                    print(f"    Alignment Sets: {alignment_sets}")
                    print(f"    Set Matches: {set_matches}")

                    matched_edges_1.remove(e1)
                    matched_edges_2.remove(e2)
                    frontier_edges_1.add(e1)
                    frontier_edges_2.add(e2)
        print(f"BEST alignment before frontier={best_alignment}")

    frontier_edges_1 = set(edges_1)
    frontier_edges_2 = set(edges_2)

    print(f"BEST alignment before back={best_alignment}\n BEST set_matches before back={best_set_matches} ")
    backtrack_vf2([], frontier_edges_1, frontier_edges_2, set(), set())
    print(f"BEST alignment after back={best_alignment}\n BEST set_matches after back={best_set_matches}")

    best_alignment, best_alignment_sets, best_set_matches = refine_alignment_with_unaligned_edges(best_alignment,
                                                                                                  best_alignment_sets,
                                                                                                  best_set_matches,
                                                                                                  sets_1, sets_2,
                                                                                                  edges_1, edges_2,
                                                                                                  current_depth)
    # set_matches, best_alignment, alignment_sets = check_implied_edge_matches(set_matches, best_alignment, alignment_sets)
    #
    #best_set_matches = finalize_isolated_set_matches(isolated_edges, best_set_matches)
    #
    # print(f"best_alignment={best_alignment}")
    # print(f"its sets={alignment_sets}")
    # print(f"its set-matches={set_matches}")
    print(f"best_alignment_after NONE={best_alignment}")
    print(f"its sets={best_alignment_sets}")
    print(f"its cost={max_cost}")
    print(f"its set-matches={best_set_matches}")

    return best_alignment, best_alignment_sets, max_cost, best_set_matches





def align_multiple_e_graphs(edges_list, sets_list):
    num_graphs = len(edges_list)
    if num_graphs < 2:
        return edges_list[0], sets_list[0], {}

    composite_edges = edges_list[0]
    composite_sets = sets_list[0]
    composite_set_matches = {}
    current_depth = 1
    for i in range(1, num_graphs):
        current_depth += 1
        print(f"[DEBUG] Aligning e-graph number {current_depth}")
        composite_edges, composite_sets, _, composite_set_matches = find_best_alignment(
            composite_edges, composite_sets, edges_list[i], sets_list[i], current_depth
        )
    return composite_edges, composite_sets, composite_set_matches




edges_list = [['e1', 'e2', 'e3', 'e4'],
                ['g1', 'g2', 'g3', 'g4'],
                ['f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7'],
['k1', 'k2', 'k3', 'k4']
                          ]
sets_list = [
    [{'e1', 'e2'}, {'e2', 'e3'}, {'e4', 'e3'}, {'e1', 'e4'}],
    [{'g1'}, {'g2', 'g3'}, {'g4', 'g3'}, {'g1', 'g2'}, {'g4'}],
[{'f1', 'f2','f5'}, {'f2', 'f3'}, {'f3', 'f4'}, {'f4','f1'}, {'f5','f6'}, {'f6','f7'}, {'f7'}],
    [{'k1'}, {'k2', 'k3'}, {'k4', 'k3'}, {'k1', 'k2'}, {'k4'}]
]




alignment_edges, alignment_sets, set_matches = align_multiple_e_graphs(edges_list, sets_list)
print("\nFinal Alignment Edges:", alignment_edges)
print("\nFinal Alignment Sets:")
for idx, edge_set in enumerate(alignment_sets):
    print(f"Set {idx + 1}: {edge_set}")

print("\nFinal Set Matches:")
for s3, (s1, s2) in set_matches.items():
    print(f" {s3}  ->  {s1},  {s2}")



