from pikachu.chem.structure import Structure
from pikachu.chem.atom import Atom
from pikachu.math_functions import Vector


class MolFileReader:
    value_to_charge = {7: -3,
                       6: -2,
                       5: -1,
                       0: 0,
                       3: 1,
                       2: 2,
                       1: 3}

    value_to_bond = {1: 'single',
                     2: 'double',
                     3: 'triple',
                     4: 'aromatic'}

    value_to_chiral_symbol = {0: None,
                              1: '/',
                              6: '\\'}

    def __init__(self, molfile=None, molfile_str=None):
        # Instantiate MolFileReader with mol_file_str or mol_file_path
        self.molfile_str = molfile_str
        self.molfile_path = molfile
        # Raise error when not enough arguments were given
        err = "molfile_str or molfile_path needed to instantiate MolFileReader"
        if not molfile and not molfile_str:
            raise ValueError(err)
        # Raise error when too many arguments were given
        if not molfile and not molfile_str:
            raise ValueError(f'{err} (both were given)')
        self.molfile_lines = self.get_molfile_lines()
        self.structure = Structure()

    def molfile_to_structure(self):
        atoms, bonds = self.get_molfile_components()
        self.parse_atom_info(atoms)
        self.parse_bond_info(bonds)
        self.structure.refine_structure()
        self.structure.set_double_bond_chirality()
        return self.structure

    def parse_atom_info(self, atoms):
        for i, atom_info in enumerate(atoms):
            atom_x = float(atom_info[:10])
            atom_y = float(atom_info[10:20])
            # atom_z = float(atom_info[20:30])
            atom_type = atom_info[31:34].strip()
            # atom_isotope = int(atom_info[34:36])
            atom_charge = self.value_to_charge[int(atom_info[36:39])]
            atom_id = i

            atom = Atom(atom_type, atom_id, None, atom_charge, False)
            atom.draw.set_position(Vector(atom_x, atom_y))
            self.structure.add_disconnected_atom(atom)

    def parse_bond_info(self, bonds):
        for i, bond_info in enumerate(bonds):
            atom_1_nr = int(bond_info[:3]) - 1
            atom_1 = self.structure.atoms[atom_1_nr]

            atom_2_nr = int(bond_info[3:6]) - 1
            atom_2 = self.structure.atoms[atom_2_nr]

            bond_type = self.value_to_bond[int(bond_info[6:9])]
            chiral_symbol = self.value_to_chiral_symbol[int(bond_info[9:12])]
            self.structure.add_bond(atom_1, atom_2, bond_type, i, chiral_symbol)

    def get_molfile_lines(self):
        # Access lines of molfile
        if self.molfile_str:
            molfile_lines = self.molfile_str.split('\n')
        else:
            with open(self.molfile_path, 'r') as mol_file:
                molfile_lines = mol_file.read()
                molfile_lines = molfile_lines.split('\n')
        return molfile_lines

    def get_molfile_components(self):
        atoms = []
        bonds = []
        # Read atom and bond counts
        line_index = 3
        counts = self.molfile_lines[line_index]
        atom_counts = int(counts[:3])
        bond_counts = int(counts[3:6])
        # Read atom information
        for _ in range(atom_counts):
            line_index += 1
            atom_info = self.molfile_lines[line_index]
            atoms.append(atom_info)
        # Read bond information
        for _ in range(bond_counts):
            line_index += 1
            bond_info = self.molfile_lines[line_index]
            bonds.append(bond_info)
        return atoms, bonds
