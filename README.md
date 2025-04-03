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
    [
      {"id": "a1", "type": "DOUBLE"}, {"id": "a2", "type": "SINGLE"}, {"id": "a3", "type": "SINGLE"},
      {"id": "a4", "type": "SINGLE"}, {"id": "a5", "type": "DOUBLE"}, {"id": "a6", "type": "SINGLE"},
      {"id": "a7", "type": "SINGLE"}
    ],
    [
      {"id": "b1", "type": "DOUBLE"}, {"id": "b2", "type": "SINGLE"}, {"id": "b3", "type": "SINGLE"},
      {"id": "b4", "type": "SINGLE"}, {"id": "b5", "type": "SINGLE"}, {"id": "b6", "type": "DOUBLE"},
      {"id": "b7", "type": "SINGLE"}, {"id": "b8", "type": "SINGLE"}, {"id": "b9", "type": "DOUBLE"},
      {"id": "b10", "type": "SINGLE"}
    ],
    [
      {"id": "c1", "type": "SINGLE"}, {"id": "c2", "type": "DOUBLE"}, {"id": "c3", "type": "SINGLE"},
      {"id": "c4", "type": "SINGLE"}, {"id": "c5", "type": "SINGLE"}, {"id": "c6", "type": "DOUBLE"}
    ],
    [
      {"id": "d1", "type": "SINGLE"}, {"id": "d2", "type": "SINGLE"}, {"id": "d3", "type": "DOUBLE"},
      {"id": "d4", "type": "SINGLE"}, {"id": "d5", "type": "SINGLE"}, {"id": "d6", "type": "SINGLE"},
      {"id": "d7", "type": "DOUBLE"}
    ],
    [
      {"id": "e1", "type": "SINGLE"}, {"id": "e2", "type": "DOUBLE"}, {"id": "e3", "type": "SINGLE"},
      {"id": "e4", "type": "SINGLE"}, {"id": "e5", "type": "SINGLE"}, {"id": "e6", "type": "DOUBLE"}
    ]
  ],
  "sets_list": [
    [
      {"edges": ["a1"], "set_type": "O"}, {"edges": ["a1", "a2", "a7"], "set_type": "C"},
      {"edges": ["a2", "a3"], "set_type": "C"}, {"edges": ["a3", "a4"], "set_type": "C"},
      {"edges": ["a4", "a5", "a6"], "set_type": "C"}, {"edges": ["a5"], "set_type": "O"},
      {"edges": ["a6", "a7"], "set_type": "N"}
    ],
    [
      {"edges": ["b1"], "set_type": "O"}, {"edges": ["b1", "b2", "b3"], "set_type": "C"},
      {"edges": ["b2"], "set_type": "O"}, {"edges": ["b3", "b4", "b10"], "set_type": "C"},
      {"edges": ["b4", "b5"], "set_type": "N"}, {"edges": ["b5", "b6", "b7"], "set_type": "C"},
      {"edges": ["b6"], "set_type": "O"}, {"edges": ["b7", "b8"], "set_type": "N"},
      {"edges": ["b8", "b9", "b10"], "set_type": "C"}, {"edges": ["b9"], "set_type": "O"}
    ],
    [
      {"edges": ["c1"], "set_type": "C"}, {"edges": ["c1", "c2", "c3"], "set_type": "C"},
      {"edges": ["c2"], "set_type": "O"}, {"edges": ["c3", "c4"], "set_type": "N"},
      {"edges": ["c4", "c5", "c6"], "set_type": "C"}, {"edges": ["c5"], "set_type": "N"},
      {"edges": ["c6"], "set_type": "O"}
    ],
    [
      {"edges": ["d1"], "set_type": "C"}, {"edges": ["d1", "d2"], "set_type": "C"},
      {"edges": ["d2", "d3", "d4"], "set_type": "C"}, {"edges": ["d3"], "set_type": "O"},
      {"edges": ["d4", "d5"], "set_type": "N"}, {"edges": ["d5", "d6", "d7"], "set_type": "C"},
      {"edges": ["d6"], "set_type": "N"}, {"edges": ["d7"], "set_type": "O"}
    ],
    [
      {"edges": ["e1"], "set_type": "N"}, {"edges": ["e1", "e2", "e3"], "set_type": "C"},
      {"edges": ["e2"], "set_type": "O"}, {"edges": ["e3", "e4"], "set_type": "N"},
      {"edges": ["e4", "e5", "e6"], "set_type": "C"}, {"edges": ["e5"], "set_type": "N"},
      {"edges": ["e6"], "set_type": "O"}
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


**Converting SMILES Strings to Input Format:**

If you're working with chemical structures, you can generate valid input JSON for `EdgeAlign` using the provided script: `smilesToInput.py`.
This script takes one or more SMILES strings and outputs a properly formatted JSON file with the required `edges_list` and `sets_list`.

```bash
python smilesToInput.py "O=C1CCC(=O)N1" "O=C1NC(=O)C(C(=O)O)N1" > input.json
```

Each SMILES string is treated as a separate molecule.
The resulting JSON can be used directly as input to `EdgeAlign`.

To generate and view a labeled image of the molecules:
python smilesToInput.py "O=C1CCC(=O)N1" "O=C1NC(=O)C(C(=O)O)N1" --win > input.json

This will:
- Save a PNG file `all_labeled_molecules.png` showing each molecule with labeled bond IDs.
- Open the image in a pop-up window (if a graphical interface is available).

📤 Output

`stdout`: JSON structure compatible with EdgeAlign
`all_labeled_molecules.png`: Visual representation of molecules with bond labels
No extra terminal output unless `--win` is used

💡 Notes
- Invalid SMILES strings are silently ignored.
- Edge labels are prefixed (e.g., `a1`, `b2`, ...) for uniqueness and clarity.
- This script is especially useful for cheminformatics workflows.


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
