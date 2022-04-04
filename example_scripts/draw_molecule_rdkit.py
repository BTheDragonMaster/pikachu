from sys import argv

from rdkit import Chem
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)


def draw_molecule(smiles, out_file):
    mol = MolFromSmiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.drawOptions().addStereoAnnotation = False
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    with open(out_file, 'w') as out:
        out.write(svg)

if __name__ == "__main__":
    smiles = argv[1]
    out_file = argv[2]
    draw_molecule(smiles, out_file)


