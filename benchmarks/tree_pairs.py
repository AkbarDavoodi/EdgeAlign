"""Appendix B benchmark: pairwise alignment of synthetic tree pairs.

Each instance pairs the path on vertices 0, ..., n with the tree obtained by
replacing the edge {n-1, n} with {1, n}.  All edge labels agree.  Vertex
labels are either uniform or repeat C/N/O in vertex-index order.  For each n
the median of three runs is reported; a series stops at its first timeout.

Usage (from the repository root):

    python benchmarks/tree_pairs.py                       # current EdgeAlign.py
    python benchmarks/tree_pairs.py --module path/to/EdgeAlign_direct.py

The direct recursive prototype compared in the paper is the EdgeAlign.py of
git revision 7b0ef3beb0ae1bddbbb1db7a0a7739a4fed273bb; extract it
from a Git clone with:
    git show 7b0ef3beb0ae1bddbbb1db7a0a7739a4fed273bb:EdgeAlign.py > EdgeAlign_direct.py
"""
import argparse
import importlib.util
import os
import signal
import statistics
import sys
import time


def tree_pair(n, labelled):
    label = (lambda v: "CNO"[v % 3]) if labelled else (lambda v: "X")

    def build(prefix, bonds):
        edges = [{"id": f"{prefix}{i}", "type": "SINGLE"} for i in range(len(bonds))]
        incident = {}
        for i, (u, v) in enumerate(bonds):
            incident.setdefault(u, []).append(edges[i]["id"])
            incident.setdefault(v, []).append(edges[i]["id"])
        sets = [{"edges": incident[v], "set_type": label(v)} for v in sorted(incident)]
        return edges, sets

    path = [(i, i + 1) for i in range(n)]
    tree = path[:-1] + [(1, n)]
    return build("p", path), build("t", tree)


class Timeout(Exception):
    pass


def load(path):
    spec = importlib.util.spec_from_file_location("edgealign_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.configure_logging("CRITICAL")
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--module", default=os.path.join(os.path.dirname(__file__), "..", "EdgeAlign.py"))
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=60, help="seconds per run")
    parser.add_argument("--max-n", type=int, default=20)
    args = parser.parse_args()
    ea = load(args.module)

    def on_alarm(*_):
        raise Timeout

    signal.signal(signal.SIGALRM, on_alarm)
    print(f"# module={args.module} python={sys.version.split()[0]} runs={args.runs} timeout={args.timeout}s")
    print("family\tn\tmedian_s\truns_s")
    for labelled in (False, True):
        family = "CNO" if labelled else "uniform"
        for n in range(3, args.max_n + 1):
            (e1, s1), (e2, s2) = tree_pair(n, labelled)
            times = []
            for _ in range(args.runs):
                signal.alarm(args.timeout)
                try:
                    start = time.perf_counter()
                    ea.align_multiple_e_graphs([e1, e2], [s1, s2])
                    times.append(time.perf_counter() - start)
                except Timeout:
                    times = None
                    break
                finally:
                    signal.alarm(0)
            if times is None:
                print(f"{family}\t{n}\ttimeout")
                break
            print(f"{family}\t{n}\t{statistics.median(times):.9f}\t"f"{' '.join(f'{t:.9f}' for t in times)}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
