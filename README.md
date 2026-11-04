# EdgeAlign

**EdgeAlign** is a Python reference implementation of the E-graph alignment framework presented in *Edge-wise Graph Alignments*. It treats graph edges as the objects being aligned, with their endpoint relationships represented by incidence sets. For molecular graphs, these correspond to bonds and atoms. The resulting alignment is itself an E-graph.

The implementation follows Algorithms 5 and 6 of the paper. Each pairwise step uses exact search for the prototype’s objective on fixed representatives. The complete multiple alignment is progressive: its result can depend on input order, and global optimality is not guaranteed.

## Quick start

Python 3.12 or newer is recommended. The core implementation requires NetworkX and Matplotlib.

```bash
git clone https://github.com/AkbarDavoodi/EdgeAlign.git
cd EdgeAlign
python -m venv .venv
source .venv/bin/activate
python -m pip install networkx matplotlib
python EdgeAlign.py --input examples/showcase_input.json --log-level INFO
```

This runs the five-molecule example used in the paper. Results are printed through logging. Add `--log-file showcase.log` to save them. The supplied `examples/showcase_expected.json` records the expected alignment and progressive scores.

## Your own inputs

Replace the example path with your JSON input file. Inputs contain `edges_list` and `sets_list`, describing labelled edges and their endpoint incidence sets. Edge identifiers must be unique across the complete input. Matches require equal edge and vertex labels.

Use [`examples/showcase_input.json`](examples/showcase_input.json) as a complete example. The [usage guide](docs/usage.md) explains the input rules, Python interface and output representation.

Optional SMILES conversion and molecular drawing helpers additionally require RDKit and Pillow; their commands and options are documented in the usage guide.

## Reproduce Figure B1

Run both implementations compared in Appendix B with:

```bash
python benchmarks/tree_pairs.py
```

The driver automatically selects the supplied **Direct recursive search** and **Algorithm 6** implementations. It reports medians of three runs with a 60-second limit per alignment call. Results are written to `tree-pairs-rerun/`; use `--output another-directory` for another run.

Both measured code versions and the original timing records are included. The command works from a clone or ZIP download. Absolute runtimes and timeout positions depend on the execution environment.

See the [benchmark documentation](benchmarks/appendix_b/README.md) for the protocol, recorded environment and original results. The benchmark runs on macOS and Linux; Windows users can use WSL.

## Tests

Run the regression checks with:

```bash
python -m unittest discover -s tests -v
```

## Documentation

The [usage guide](docs/usage.md) covers detailed input requirements, command-line options, the Python API, molecular conversion and drawing tools, and the implementation’s objective and limitations.

## Citation and license

Please cite:

> Akbar Davoodi, Christoph Flamm, Daniel Merkle, and Peter F. Stadler. *Edge-wise Graph Alignments*. Submitted.

EdgeAlign is distributed under the [MIT License](LICENSE).
