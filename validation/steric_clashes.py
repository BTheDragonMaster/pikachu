from sys import argv

from pikachu.general import read_smiles, position_smiles
from pikachu.drawing.drawing import Drawer
from pikachu.math_functions import Vector

from rdkit.Chem.rdCoordGen import AddCoords
from rdkit.Chem import MolFromSmiles

import timeout_decorator, time


def get_rdkit_coords(smiles):
    structure = MolFromSmiles(smiles)
    AddCoords(structure)
    conformer = structure.GetConformer()
    atoms = structure.GetAtoms()

    atom_indices = []

    atom_positions = []
    for atom in atoms:
        position = conformer.GetAtomPosition(atom.GetIdx())
        atom_positions.append(Vector(position.x, position.y))
        atom_indices.append(atom.GetIdx())

    bond_lookup = {}

    for bond in structure.GetBonds():
        atom_1 = bond.GetBeginAtomIdx()
        atom_2 = bond.GetEndAtomIdx()
        if atom_1 not in bond_lookup:
            bond_lookup[atom_1] = {}
        if atom_2 not in bond_lookup:
            bond_lookup[atom_2] = {}

        bond_lookup[atom_1][atom_2] = bond
        bond_lookup[atom_2][atom_1] = bond

    return atom_indices, atom_positions, bond_lookup


def get_pikachu_coords(smiles):
    structure = read_smiles(smiles)
    drawing = position_smiles(smiles)
    atom_positions = []
    atoms = []
    for atom in drawing.structure.graph:
        if atom.draw.is_drawn:
            atoms.append(atom)
            atom_positions.append(atom.draw.position)

    return atoms, atom_positions, structure.bond_lookup


def find_average_bond_length(atoms, atom_positions, bond_lookup):
    counter = 0
    total_bond_length = 0

    for i, atom_1 in enumerate(atoms):
        atom_position_1 = atom_positions[i]
        for j, atom_2 in enumerate(atoms):
            atom_position_2 = atom_positions[j]
            if i != j:
                if atom_1 in bond_lookup and atom_2 in bond_lookup[atom_1]:
                    bond_length = atom_position_1.get_distance(atom_position_2)
                    total_bond_length += bond_length
                    counter += 1

    try:
        average_bond_length = total_bond_length / float(counter)
    except ZeroDivisionError:
        return None

    return average_bond_length


def find_steric_clashes(atoms, atom_positions, average_bond_length, bond_lookup):
    steric_clashes = 0
    for i, position_1 in enumerate(atom_positions):
        atom_1 = atoms[i]
        for j, position_2 in enumerate(atom_positions):
            atom_2 = atoms[j]
            if atom_2 in bond_lookup:
                if i != j and atom_1 not in bond_lookup[atom_2]:
                    distance = position_1.get_distance(position_2)
                    if distance < 0.5 * average_bond_length:
                        steric_clashes += 1
            elif i != j:
                distance = position_1.get_distance(position_2)
                if distance < 0.5 * average_bond_length:
                    steric_clashes += 1

    return steric_clashes


# @timeout_decorator.timeout(40)
def find_clashes_rdkit(smiles):
    atoms, atom_positions, bond_lookup = get_rdkit_coords(smiles)
    av_bond_length = find_average_bond_length(atoms, atom_positions, bond_lookup)

    if av_bond_length:
        clashes = find_steric_clashes(atoms, atom_positions, av_bond_length, bond_lookup)


    # Only happens if there are no bonds in the structure.
    else:
        clashes = 0
    if clashes:
        print(f"Handling smiles: {smiles}, found {clashes} clashes (RDKit).")
    return clashes


# @timeout_decorator.timeout(40)
def find_clashes_pikachu(smiles):
    atoms, atom_positions, bond_lookup = get_pikachu_coords(smiles)
    av_bond_length = find_average_bond_length(atoms, atom_positions, bond_lookup)

    if av_bond_length:
        clashes = find_steric_clashes(atoms, atom_positions, av_bond_length, bond_lookup)

    # Only happens if there are no bonds in the structure.
    else:
        clashes = 0

    if clashes:
        print(f"Handling smiles: {smiles}, found {clashes} clashes (PIKAChU).")
    return clashes


def assess_tools(smiles_file, failed_out, clashes_out):
    steric_clashes_rdkit = 0
    steric_clashes_pikachu = 0

    clashing_structures_rdkit = 0
    clashing_structures_pikachu = 0

    failed_smiles_rdkit = 0
    failed_smiles_pikachu = 0
    clashes = open(clashes_out, 'w')
    failed = open(failed_out, 'w')
    with open(smiles_file, 'r') as smi:
        for i, line in enumerate(smi):
            smiles = line.strip()
            # try:
            #     clashes_rdkit = find_clashes_rdkit(smiles)
            #     steric_clashes_rdkit += clashes_rdkit
            #
            #     if clashes_rdkit:
            #         clashing_structures_rdkit += 1
            #
            # except Exception:
            #     print(f"Failed smiles RDKIT: {smiles}")
            #     failed_smiles_rdkit += 1

            try:
                clashes_pikachu = find_clashes_pikachu(smiles)
                steric_clashes_pikachu += clashes_pikachu

                if clashes_pikachu:
                    clashes.write(f"{smiles}\t{clashes_pikachu}\n")
                    clashing_structures_pikachu += 1

            except Exception as e:
                print(f"Failed smiles PIKAChU: {smiles}")
                failed_smiles_pikachu += 1
                failed.write(f"{smiles}\n")

            if i % 500 == 0:
                print(f"Processed {i} smiles.")
                # print(f"Steric clashes RDKit: {steric_clashes_rdkit / 2}")
                print(f"Steric clashes PIKAChU: {steric_clashes_pikachu / 2}")
                # print(f"Clashing structures RDKit: {clashing_structures_rdkit}")
                print(f"Clashing structures PIKAChU: {clashing_structures_pikachu}")
                # print(f"Failed SMILES RDKit: {failed_smiles_rdkit}")
                print(f"Failed SMILES PIKAChU: {failed_smiles_pikachu}")
                print('\n')
                
            if i == 100000:

                break
    clashes.close()
    failed.close()

    # print(f"Steric clashes RDKit: {steric_clashes_rdkit / 2}")
    print(f"Steric clashes PIKAChU: {steric_clashes_pikachu / 2}")
    # print(f"Clashing structures RDKit: {clashing_structures_rdkit}")
    print(f"Clashing structures PIKAChU: {clashing_structures_pikachu}")
    # print(f"Failed SMILES RDKit: {failed_smiles_rdkit}")
    print(f"Failed SMILES PIKAChU: {failed_smiles_pikachu}")


if __name__ == "__main__":
    smiles_file = argv[1]
    failed = smiles_file.split('.')[0] + '_failed_steric.txt'
    clashing = smiles_file.split('.')[0] + '_clashing_steric.txt'
    assess_tools(smiles_file, failed, clashing)




