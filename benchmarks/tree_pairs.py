"""Appendix B tree pairs, with the original inputs and timing protocol.

Run from the repository root. See benchmarks/appendix_b/README.md for the
measured revisions, original records, and commands for comparing both versions.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import signal
import statistics
import subprocess
import sys
import time

from check_alignment import validate_output


def tree_pair(n, labelled):
    graphs = []
    for row in range(2):
        bonds = [(v, v + 1) for v in range(n)]
        if row:
            bonds[-1] = (1, n)
        edges = [{"id": f"g{row}e{i}", "type": "SINGLE"} for i in range(n)]
        sets = [{"edges": [edges[i]["id"] for i, bond in enumerate(bonds) if v in bond],
                 "set_type": "CNO"[v % 3] if labelled else "C"} for v in range(n + 1)]
        graphs.append((edges, sets))
    return graphs


def task(n, labelled):
    return dict(id=f'pair_{"CNO" if labelled else "C"}_{n}', n=n,
                family="CNO" if labelled else "uniform", expected=n - 1,
                inputs=list(map(list, zip(*tree_pair(n, labelled)))))


def worker(module_path, n, labelled, limit):
    spec = importlib.util.spec_from_file_location("edgealign_under_test", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.configure_logging("CRITICAL")
    case = task(n, labelled)

    def on_alarm(*_):
        raise TimeoutError("alignment-call limit")

    signal.signal(signal.SIGALRM, on_alarm)
    started = time.perf_counter()
    signal.setitimer(signal.ITIMER_REAL, limit)
    try:
        output = module.align_multiple_e_graphs(*case["inputs"])
        seconds = time.perf_counter() - started
    except TimeoutError:
        return dict(status="timeout", seconds=time.perf_counter() - started)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    check = validate_output(output, case["inputs"])
    if not check["valid"] or check["core_bonds"] != case["expected"]:
        raise ValueError(f"Invalid or suboptimal benchmark result: {check}")
    return dict(status="complete", seconds=seconds, valid=True, score=check["core_bonds"])


def summarize(rows):
    summary = []
    for method, case in dict.fromkeys((r["method"], r["case"]) for r in rows):
        own = [r for r in rows if (r["method"], r["case"]) == (method, case)]
        complete = [r["seconds"] for r in own if r["status"] == "complete"]
        timeout = any(r["status"] == "timeout" for r in own)
        summary.append(dict(method=method, case=case, n=own[0]["n"], family=own[0]["family"],
                            completed=len(complete), attempts=len(own),
                            median_seconds=None if timeout else statistics.median(complete),
                            timeout=timeout))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", type=Path, default=Path(__file__).resolve().parents[1] / "EdgeAlign.py")
    parser.add_argument("--baseline", type=Path, help="direct recursive version; enables the matched comparison")
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=60, help="seconds per alignment call")
    parser.add_argument("--max-n", type=int, default=40)
    parser.add_argument("--output", type=Path, default=Path("tree-pairs-rerun"), help="new results directory; never overwritten")
    args = parser.parse_args()
    if args.runs < 1 or args.timeout <= 0 or args.max_n < 3:
        parser.error("require --runs >= 1, --timeout > 0, and --max-n >= 3")
    modules = {"algorithm6": args.module.resolve()}
    if args.baseline:
        modules = {"before": args.baseline.resolve(), **modules}
    for path in modules.values():
        if not path.is_file():
            parser.error(f"module does not exist: {path}")
    args.output.mkdir(parents=True, exist_ok=False)

    def save(name, value):
        (args.output / name).write_text(json.dumps(value, indent=2) + "\n")

    save("protocol.json", dict(limit=args.timeout, repeats=args.runs, python=sys.version,
         platform=platform.platform(), machine=platform.machine(), hash_seed=0,
         timing="alignment call only; imports and independent validation excluded",
         order="serial; before/algorithm6 alternate per repeat through n=8; extend algorithm6, then before",
         stopping="stop repetitions and the method/family at its first timeout",
         sources={name: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                  for name, path in modules.items()}))
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg")
    rows, tasks, stopped = [], {}, set()
    with (args.output / "runs.jsonl").open("w") as stream:
        def measure(method, n, labelled, repeat):
            case = task(n, labelled)
            tasks[case["id"]] = case
            call = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker",
                                   str(modules[method]), str(n), str(int(labelled)), str(args.timeout)],
                                  capture_output=True, text=True, env=env, timeout=args.timeout + 20)
            if call.returncode:
                raise RuntimeError(call.stderr)
            row = json.loads(call.stdout)
            row.update(method=method, case=case["id"], family=case["family"], n=n, repeat=repeat)
            rows.append(row)
            stream.write(json.dumps(row) + "\n")
            stream.flush()
            save("inputs.json", list(tasks.values()))
            save("summary.json", summarize(rows))
            print(method, case["id"], repeat, row["status"], f'{row["seconds"]:.6f}s', flush=True)
            if row["status"] == "timeout":
                stopped.add((method, labelled))

        for n in range(3, min(8, args.max_n) + 1):
            for labelled in (False, True):
                for repeat in range(args.runs):
                    order = list(modules) if repeat % 2 == 0 else list(reversed(modules))
                    for method in order:
                        if (method, labelled) not in stopped:
                            measure(method, n, labelled, repeat)
        for method in ("algorithm6", "before"):
            if method not in modules:
                continue
            for labelled in (False, True):
                for n in range(9, args.max_n + 1):
                    if (method, labelled) in stopped:
                        break
                    for repeat in range(args.runs):
                        if (method, labelled) not in stopped:
                            measure(method, n, labelled, repeat)
    print("method\tfamily\tn\tmedian_s")
    for row in summarize(rows):
        value = "timeout" if row["timeout"] else f'{row["median_seconds"]:.9f}'
        print(f'{row["method"]}\t{row["family"]}\t{row["n"]}\t{value}')


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        print(json.dumps(worker(sys.argv[2], int(sys.argv[3]), bool(int(sys.argv[4])), float(sys.argv[5]))))
    else:
        main()
