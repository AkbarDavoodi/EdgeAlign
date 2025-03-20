import sys
import json
import io
from PIL import Image
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

def int_to_alpha(num):
    """
    Convert a 1-based integer to a lexicographic string:
       1 -> 'a'
       2 -> 'b'
       ...
       26 -> 'z'
       27 -> 'aa'
       28 -> 'ab'
       etc.
    """
    s = []
    while num > 0:
        num -= 1
        s.append(chr(ord('a') + (num % 26)))
        num //= 26
    return ''.join(reversed(s))

def generate_graph_lists(mol, mol_index):
    """
    Given an RDKit Mol, returns two lists:
     - edges_list: [{"id": bond_id, "type": bond_type}, ...]
     - sets_list: [{"edges": [...], "set_type": <atom symbol>}, ...]

    Uses a lexicographic prefix for bond IDs (a, b, ... z, aa, ab, ...),
    based on the 1-based mol_index.
    """
    rdDepictor.Compute2DCoords(mol)
    prefix = int_to_alpha(mol_index)

    bond_counter = 1
    edges_list = []
    sets_list = []
    edge_id_map = {}

    # Build edges_list
    for bond in mol.GetBonds():
        bond_id = f"{prefix}{bond_counter}"
        bond_counter += 1
        
        bond_type = str(bond.GetBondType())  # e.g. 'SINGLE', 'DOUBLE'
        edges_list.append({"id": bond_id, "type": bond_type})

        # Auto-label in the RDKit drawing
        bond.SetProp("bondNote", bond_id)

        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        key = tuple(sorted((a1, a2)))
        edge_id_map[key] = bond_id

    # Build sets_list (one entry per atom)
    for atom in mol.GetAtoms():
        incident_edges = []
        for bond in atom.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            key = tuple(sorted((a1, a2)))
            if key in edge_id_map:
                incident_edges.append(edge_id_map[key])
        sets_list.append({
            "edges": incident_edges,
            "set_type": atom.GetSymbol()
        })

    return edges_list, sets_list


def draw_molecule_as_pil(mol, size=(400, 400), annotation_scale=0.5):
    """
    Draw the RDKit Mol into a PIL image (PNG) with bond labels (bondNote).
    Returns the PIL image.
    """
    width, height = size
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.annotationFontScale = annotation_scale
    
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    png_data = drawer.GetDrawingText()
    return Image.open(io.BytesIO(png_data))


def combine_images_side_by_side(images, padding=10, background_color=(255, 255, 255)):
    """
    Places a list of PIL images side by side into one wide image.
    Returns the combined PIL image, or None if no images provided.
    """
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


def main():
    args = sys.argv[1:]

    # Detect optional --win flag
    show_window = False
    if "--win" in args:
        show_window = True
        args.remove("--win")

    # If no SMILES after removing --win, just output empty JSON
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
            # Skip invalid SMILES silently
            continue

        edges_list, sets_list = generate_graph_lists(mol, mol_index)
        all_edges.append(edges_list)
        all_sets.append(sets_list)

        # Draw the molecule into a PIL image
        pil_img = draw_molecule_as_pil(mol, size=(400, 400), annotation_scale=0.5)
        molecule_images.append(pil_img)

        mol_index += 1

    # Prepare final JSON
    result = {
        "edges_list": all_edges,
        "sets_list": all_sets
    }

    # Print only the JSON, no extra messages
    print(json.dumps(result, indent=2))

    # If no valid molecules or no images, we're done
    if not molecule_images:
        return

    # Create one side-by-side collage
    collage = combine_images_side_by_side(molecule_images, padding=20)
    if collage is None:
        return

    # Save collage to a file (silent)
    collage_filename = "all_labeled_molecules.png"
    collage.save(collage_filename)

    # If --win was given, display that collage in a window
    if show_window:
        plt.figure()
        plt.imshow(collage)
        plt.axis("off")
        plt.show()


if __name__ == "__main__":
    main()


