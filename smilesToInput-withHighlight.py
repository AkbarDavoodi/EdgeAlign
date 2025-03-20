import sys
import json
import io
import ast
from PIL import Image
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


###############################################################################
# 1) Utility to flatten nested solution tuples, e.g. ((((None, 'b1'), 'c1'), 'd5'), 'e4')
###############################################################################
def flatten_nested_tuples(obj):
    """
    Recursively flattens a nested tuple structure of the form:
      ((((None, 'b1'), 'c1'), 'd5'), 'e4')
    into a Python list [None, 'b1', 'c1', 'd5', 'e4'].
    """
    if obj is None or isinstance(obj, str):
        return [obj]
    out = []
    for x in obj:
        out.extend(flatten_nested_tuples(x))
    return out

###############################################################################
# 2) Parse the solution file. We store how many graphs each edge ID is mapped in.
###############################################################################
def parse_solution_file(solution_path):
    """
    Looks for lines beginning with 'INFO:'
    Each line has a nested tuple after 'INFO:' that we parse via `ast.literal_eval`.
    We flatten the tuple, remove None, and see how many edges appear.
    If at least 2 edges appear, we record that count (2..5) in highlight_count[edge_id].
    If the same edge ID appears in multiple lines with different cardinalities, we keep
    the maximum.
    """
    highlight_count = {}
    
    with open(solution_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("INFO:"):
                continue
            nested_str = line[5:].strip()
            try:
                nested_obj = ast.literal_eval(nested_str)
            except Exception:
                continue
            flat = flatten_nested_tuples(nested_obj)
            edge_ids = [x for x in flat if x is not None]
            count = len(edge_ids)
            if count < 2:
                continue  # Only care about edges that appear in 2..5 graphs

            for e_id in edge_ids:
                prev = highlight_count.get(e_id, 0)
                if count > prev:
                    highlight_count[e_id] = count
    return highlight_count

###############################################################################
# 3) Convert integer -> lexicographic prefix, e.g. 1->'a', 2->'b', ..., 26->'z', 27->'aa'
###############################################################################
def int_to_alpha(num):
    s = []
    while num > 0:
        num -= 1
        s.append(chr(ord('a') + (num % 26)))
        num //= 26
    return ''.join(reversed(s))

###############################################################################
# 4) Generate edges_list + sets_list for an RDKit Mol. Also store bond IDs for labeling.
###############################################################################
def generate_graph_lists(mol, mol_index):
    rdDepictor.Compute2DCoords(mol)

    prefix = int_to_alpha(mol_index)
    bond_counter = 1

    edges_list = []
    sets_list = []
    edge_id_map = {}

    for bond in mol.GetBonds():
        bond_id = f"{prefix}{bond_counter}"
        bond_counter += 1

        bond_type = str(bond.GetBondType())
        edges_list.append({"id": bond_id, "type": bond_type})

        # Let RDKit auto-draw this label
        bond.SetProp("bondNote", bond_id)

        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        key = tuple(sorted((a1, a2)))
        edge_id_map[key] = bond_id

    for atom in mol.GetAtoms():
        incident_edges = []
        for b in atom.GetBonds():
            a1 = b.GetBeginAtomIdx()
            a2 = b.GetEndAtomIdx()
            key = tuple(sorted((a1, a2)))
            if key in edge_id_map:
                incident_edges.append(edge_id_map[key])
        sets_list.append({
            "edges": incident_edges,
            "set_type": atom.GetSymbol()
        })

    return edges_list, sets_list

###############################################################################
# 5) Draw the molecule as a PIL image, highlighting edges that appear in 2..5 graphs.
###############################################################################
def draw_molecule_as_pil(mol, mol_index, highlight_count, size=(400, 400), annotation_scale=0.5):
    """
    We'll highlight each bond based on how many graphs it appears in:
      5 -> red
      4 -> orange
      3 -> pink
      2 -> light grey
    Uses the older DrawMolecule() approach with highlightBonds & highlightBondColors.
    """
    # Lexicographic prefix for this molecule
    prefix = int_to_alpha(mol_index)

    color_map = {
        2: (0.83, 0.83, 0.83),  # light grey
        3: (1.0, 0.75, 0.8),    # pink
        4: (1.0, 0.65, 0.0),    # orange
        5: (1.0, 0.0, 0.0),     # red
    }

    highlight_bonds = []
    highlight_bond_colors = {}

    # Identify which bonds need highlighting
    for b_idx, bond in enumerate(mol.GetBonds()):
        bond_id = bond.GetProp("bondNote") if bond.HasProp("bondNote") else None
        if bond_id and bond_id.startswith(prefix):
            c = highlight_count.get(bond_id, 0)
            # highlight only if c >= 2
            if c >= 2:
                if c > 5:
                    c = 5
                color = color_map[c]
                highlight_bonds.append(b_idx)
                highlight_bond_colors[b_idx] = color

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.annotationFontScale = annotation_scale
    
    # 'highlightBonds' is a list of bond indices to highlight
    # 'highlightBondColors' is a dict {bond_idx: (r,g,b)}
    drawer.DrawMolecule(
        mol,
        highlightAtoms=[],
        highlightBonds=highlight_bonds,
        highlightAtomColors={},
        highlightBondColors=highlight_bond_colors
    )
    drawer.FinishDrawing()

    png_data = drawer.GetDrawingText()
    return Image.open(io.BytesIO(png_data))

###############################################################################
# 6) Combine images side by side into a single PNG
###############################################################################
def combine_images_side_by_side(images, padding=10, background_color=(255, 255, 255)):
    if not images:
        return None

    widths = [img.width for img in images]
    heights = [img.height for img in images]

    total_width = sum(widths) + padding * (len(images) - 1)
    max_height = max(heights)

    combined_img = Image.new('RGB', (total_width, max_height), background_color)

    x_offset = 0
    for img in images:
        combined_img.paste(img, (x_offset, 0))
        x_offset += img.width + padding

    return combined_img

###############################################################################
# 7) Main script logic
###############################################################################
def main():
    args = sys.argv[1:]
    show_window = False
    highlight_count = {}

    # Check if --win is specified
    if "--win" in args:
        show_window = True
        args.remove("--win")

    # Check for --sol=filename or --sol filename
    solution_file = None
    for i, arg in enumerate(args):
        if arg is None:
            continue
        if arg.startswith("--sol="):
            solution_file = arg.split("=", 1)[1]
            args[i] = None
        elif arg == "--sol" and i + 1 < len(args):
            solution_file = args[i + 1]
            args[i] = None
            args[i + 1] = None
            break

    # Clean up None's from arg list
    args = [a for a in args if a is not None]

    # If we have a solution file, parse it
    if solution_file:
        highlight_count = parse_solution_file(solution_file)

    # If no SMILES are given after removing flags, print empty JSON
    if not args:
        result = {"edges_list": [], "sets_list": []}
        print(json.dumps(result))
        return

    all_edges = []
    all_sets = []
    molecule_images = []

    mol_index = 1
    for smi in args:
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            # skip invalid SMILES
            continue

        edges_list, sets_list = generate_graph_lists(mol, mol_index)
        all_edges.append(edges_list)
        all_sets.append(sets_list)

        # Draw image with highlights
        img = draw_molecule_as_pil(
            mol,
            mol_index,
            highlight_count,
            size=(400, 400),
            annotation_scale=0.5
        )
        molecule_images.append(img)

        mol_index += 1

    # Prepare final JSON
    result = {
        "edges_list": all_edges,
        "sets_list": all_sets
    }
    # Print JSON only
    print(json.dumps(result, indent=2))

    # If no valid molecules, we're done
    if not molecule_images:
        return

    # Combine images side by side
    collage = combine_images_side_by_side(molecule_images, padding=20)
    if collage is None:
        return

    # Save the collage (silent)
    collage_filename = "all_labeled_molecules.png"
    collage.save(collage_filename)

    # If --win is given, show the collage in a window
    if show_window:
        plt.figure()
        plt.imshow(collage)
        plt.axis("off")
        plt.show()

if __name__ == "__main__":
    main()


