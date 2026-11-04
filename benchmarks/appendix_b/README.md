# Original Appendix B measurements

These are the original observations used in the paper, copied without changes
from the archived scaling study. No new timing campaign is included.

| File | Contents |
| --- | --- |
| `inputs.json` | The 22 exact tree-pair inputs, including identifiers, labels and expected optimum. |
| `runs.jsonl` | All 100 attempts: 96 completed calls and four timeouts; elapsed seconds and repetition indices. |
| `summary.json` | The 36 method/instance summaries; medians require three completed calls. |
| `protocol.json` | Original measurement protocol, including Python version, timing scope and source hashes. |
| `versions.json` | Full Git commit IDs and SHA-256 hashes of the two measured `EdgeAlign.py` files. |
| `environment.json` | Measurement Python version and recovered machine/software details, with provenance. |

The method names `before` and `algorithm6` in these records mean **Direct
recursive search** and **Endpoint-orientation search (Alg. 6)**, respectively.
The former includes correctness fixes; it is not the unmodified submitted code.
Its measured commit is `7b0ef3beb0ae1bddbbb1db7a0a7739a4fed273bb`; the Algorithm 6
commit is `ceaa0db28806dffd5a0c09286fe9e97b526a233a`. Both files are available
in the repository history; `EdgeAlign0.py` was not used in this comparison.

## Measurement protocol

Each pair comprises a path on vertices `0, ..., n` and a second tree obtained
by replacing `{n-1, n}` with `{1, n}`. Edges have label `SINGLE`, identifiers
`g0e0`, ... and `g1e0`, ...; vertex labels are `C` or `CNO[v % 3]`.
The known optimum is `n-1` matched edges.

Every call ran in a fresh process with `PYTHONHASHSEED=0`,
`PYTHONDONTWRITEBYTECODE=1` and `MPLBACKEND=Agg`. The 60-second limit and timing
cover only `align_multiple_e_graphs`; imports, input construction and independent
result validation are excluded. Each completed instance has three repetitions.
The first timeout stops further repetitions and extension of that method/family.
Timeout records retain elapsed time but do not contribute a value of 60 seconds
to the plotted medians.

For sizes 3 through 8, methods alternated within each repetition, serially.
Algorithm 6 was then extended from size 9 until the first timeout in each family.
The direct recursive C/N/O series was extended afterwards with the same protocol.
The driver follows this order when `--baseline` is supplied. On a faster machine,
it can also extend the direct uniform series beyond size 8 if no timeout occurs.

The original protocol records **Python 3.13.7**. The earlier same-machine
environment record identifies an **Apple M4 Pro, 64 GiB RAM, macOS 26.6.2/arm64**
and the package versions listed in `environment.json` (including NetworkX 3.7
and Matplotlib 3.11.2). Hardware, OS and packages were not independently captured
again during this timing campaign; their archival provenance is explicit.
The README's Python 3.12.14 compatibility check is a separate observation.

## Repeating the comparison

Use a full clone, Python 3.13.7 and the archived package versions to reconstruct
the recorded environment. From the repository root:

```bash
benchmark_dir=$(mktemp -d)
git show 7b0ef3beb0ae1bddbbb1db7a0a7739a4fed273bb:EdgeAlign.py > "$benchmark_dir/EdgeAlign_direct.py"
git show ceaa0db28806dffd5a0c09286fe9e97b526a233a:EdgeAlign.py > "$benchmark_dir/EdgeAlign_algorithm6.py"
PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg python benchmarks/tree_pairs.py \
  --baseline "$benchmark_dir/EdgeAlign_direct.py" --module "$benchmark_dir/EdgeAlign_algorithm6.py" \
  --output "$benchmark_dir/results"
```

The output directory must not exist. It receives fresh `runs.jsonl`,
`inputs.json`, `summary.json` and `protocol.json`; the archived observations
remain untouched. The driver prints the medians and timeout status. Every
completed result is checked outside the timer using the original independent
projection/label checks and the known optimum, without requiring RDKit.
Absolute timings and timeout boundaries can vary with the execution environment.
