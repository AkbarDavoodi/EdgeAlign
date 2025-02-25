# EdgeAlign

**EdgeAlign** is a Python tool for aligning multiple E-graphs (edge-centric graph representations) based on the framework presented in:

> Akbar Davoodi, Christoph Flamm, Daniel Merkle, and Peter F Stadler, "Edge‑wise Graph Alignments", *Submitted*.

This repository provides a working implementation of progressive E-graph alignments using a VF2-inspired backtracking approach. The method computes pairwise alignments and then combines them into a composite multiple alignment, making it particularly useful in applications such as cheminformatics.



## Overview

In traditional graph alignment, graphs are compared based on their vertices. **EdgeAlign** represents graphs as E-graphs, where edges are the primary objects and vertices are encoded implicitly via incidence sets. This representation offers several advantages:
- **Recursive Backtracking:** Inspired by the VF2 algorithm, it recursively explores candidate edge matches.
- **Progressive Alignment:** Multiple E-graphs are aligned progressively to construct a composite alignment.




## Installation

### Prerequisites

- **Python:** Version 3.7 or later.
- **Required Packages:**  
  - `networkx`
  - `matplotlib`

Install the required packages using pip (if not already installed):

```bash
pip install networkx matplotlib
```

Then, change into the repository directory:

```bash
cd EdgeAlign
```

## Usage

The main script is `EdgeAlign.py`, which supports configurable logging and JSON input.

### Running from the Command Line

If you have a JSON input file (see [Input Format](#input-format) below), run the tool as follows:

```bash
python EdgeAlign.py --input data.json --log-level INFO
```

For more detailed (debug) output, use:

```bash
python EdgeAlign.py --input data.json --log-level DEBUG
```
If no input file is provided, the tool uses built-in test data.


### Running in an IDE (e.g., PyCharm)

To run the script in PyCharm:

1. Open **Run > Edit Configurations...** in PyCharm.
2. In the **Script parameters** field, enter:
```bash
--input data.json --log-level INFO
```
3. Ensure that `data.json` is in the same directory as `EdgeAlign.py` (or adjust the path accordingly).

## Input Format

EdgeAlign expects the input to be provided as a JSON file with two keys: `edges_list` and `sets_list`.

### Example `data.json`

```json
{
  "edges_list": [
    ["e1", "e2", "e3", "e4"],
    ["g1", "g2", "g3", "g4"],
    ["f1", "f2", "f3", "f4", "f5", "f6", "f7"],
    ["k1", "k2", "k3", "k4"]
  ],
  "sets_list": [
    [
      ["e1", "e2"],
      ["e2", "e3"],
      ["e4", "e3"],
      ["e1", "e4"]
    ],
    [
      ["g1"],
      ["g2", "g3"],
      ["g4", "g3"],
      ["g1", "g2"],
      ["g4"]
    ],
    [
      ["f1", "f2", "f5"],
      ["f2", "f3"],
      ["f3", "f4"],
      ["f4", "f1"],
      ["f5", "f6"],
      ["f6", "f7"],
      ["f7"]
    ],
    [
      ["k1"],
      ["k2", "k3"],
      ["k4", "k3"],
      ["k1", "k2"],
      ["k4"]
    ]
  ]
}
```

**Important Validations:**
- **Edge Appearance:**  
  Each edge from `edges_list` must appear in either 0 or 2 incidence sets.
  If an edge appears in any number of sets other than 0 or 2, an error is raised.
- **Distinct Edge Labels (Recommendation):**  
  It is recommended that edge labels are distinct across E-graphs for clearer output and readability.  
  Duplicate edge labels are allowed and the tool will still run correctly, but using unique labels is advised to avoid ambiguous or less interpretable results.



## How to Cite

If you use **EdgeAlign** in your research, please consider citing the following work:

> Akbar Davoodi, Christoph Flamm, Daniel Merkle, and Peter F Stadler, "Edge‑wise Graph Alignments", *Submitted*.

You can also use the following BibTeX entry:

```bibtex
@article{EdgeAlign25,
  title={Edge‑wise Graph Alignments},
  author={Davoodi, Akbar and Flamm, Christoph and Merkle, Daniel and Stadler, Peter F},
  journal={Submitted},
  year={2025}
}
```

For further theoretical details and applications, please refer to the published paper.

---

## Contributing

Contributions to **EdgeAlign** are welcome! To contribute:
- **Fork the repository** and submit pull requests.
- For major changes, please open an issue first to discuss your ideas.
- Follow the code style and include tests when applicable.

## License

This project is licensed under the [MIT License](LICENSE). See the LICENSE file for more details.

