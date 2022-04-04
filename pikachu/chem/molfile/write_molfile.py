import pkg_resources
from pikachu.drawing.drawing import Drawer, Options, draw_multiple
from pikachu.math_functions import Vector
import datetime


class MolFileWriter:
    """
    NOTE: This MolFileWriter purely exports atom coordinates (x and y) found by PIKAChU; these are unitless and should
    not be interpreted as angstrom.
    """

    charge_to_value = {-3: 7,
                       -2: 6,
                       -1: 5,
                       0: 0,
                       1: 3,
                       2: 2,
                       3: 1}

    bond_to_value = {'single': 1,
                     'double': 2,
                     'triple': 3,
                     'aromatic': 4}

    chiral_symbol_to_value = {None: 0,
                              '/': 1,
                              '\\': 6}

    def __init__(self, structure, filename, drawing_options=None, multiple=False):
        self.original_structure = structure
        if not drawing_options:
            if multiple:
                self.drawing = draw_multiple(structure, coords_only=True)
            else:
                self.drawing = Drawer(structure, coords_only=True)
        else:
            if multiple:
                self.drawing = draw_multiple(structure, coords_only=True, options=drawing_options)
            else:
                self.drawing = Drawer(structure, coords_only=True, options=drawing_options)
        self.drawn_structure = self.drawing.structure
        self.filename = filename
        self.title = filename.split('.')[0]
        self.atom_to_coords = self.get_atom_coords()
        self.datetime = datetime.datetime.now()
        self.software_version = pkg_resources.get_distribution('pikachu-chem').version
        self.atom_count = self.get_atom_count()
        self.bond_count, self.drawn_bonds = self.get_bond_count()

    def get_atom_coords(self):

        atom_to_coords = {}
        for atom in self.drawn_structure.graph:
            if atom.draw.is_drawn:
                atom_to_coords[atom] = atom.draw.position

        return atom_to_coords

    def get_atom_count(self):
        count = 0
        for atom in self.drawn_structure.graph:
            original_atom = self.original_structure.atoms[atom.nr]
            if atom.draw.is_drawn:
                count += 1
            elif original_atom.type == 'H' and original_atom.has_neighbour('N') and original_atom.get_neighbour('N').pyrrole:
                count += 1

        return count

    def get_bond_count(self):
        count = 0
        bonds = []
        for bond_nr, bond in self.drawn_structure.bonds.items():
            original_bond = self.original_structure.bonds[bond_nr]
            if (bond.atom_1.draw.is_drawn and bond.atom_2.draw.is_drawn)\
                        or (original_bond.atom_1.pyrrole and original_bond.atom_2.type == 'H')\
                        or (original_bond.atom_2.pyrrole and original_bond.atom_1.type == 'H'):
                count += 1
                bonds.append(bond)

        return count, bonds

    def write_mol_file(self):
        atom_to_line_nr = {}
        with open(self.filename, 'w') as molfile:
            molfile.write(f'{self.title}\n')
            molfile.write(f'  PIKAChU {self.software_version} {self.datetime}\n')
            molfile.write('\n')
            molfile.write(f'{str(self.atom_count).rjust(3)}{str(self.bond_count).rjust(3)}  0  0  1  0  0  0  0  0999 V2000\n')
            line_nr = 0
            for atom in self.drawn_structure.graph:
                original_atom = self.original_structure.atoms[atom.nr]
                if atom.draw.is_drawn:
                    line_nr += 1
                    atom_to_line_nr[atom] = line_nr
                    x_string = f'{atom.draw.position.x:.4f}'.rjust(10)
                    y_string = f'{atom.draw.position.y:.4f}'.rjust(10)
                    z_string = f'    0.0000'
                    charge_string = f'{str(self.charge_to_value[atom.charge]).rjust(3)}'

                    molfile.write(f'{x_string}{y_string}{z_string} {atom.type.ljust(3)} 0{charge_string}  0  0  0  0  0  0  0  0  0  0\n')
                elif original_atom.type == 'H' and original_atom.has_neighbour('N') and original_atom.get_neighbour('N').pyrrole:
                    line_nr += 1
                    atom_to_line_nr[atom] = line_nr
                    position = Vector.add_vectors(atom.get_neighbour('N').draw.position, Vector(0, -15))
                    x_string = f'{position.x:.4f}'.rjust(10)
                    y_string = f'{position.y:.4f}'.rjust(10)
                    z_string = f'    0.0000'
                    molfile.write(
                        f'{x_string}{y_string}{z_string} {atom.type.ljust(3)} 0{charge_string}  0  0  0  0  0  0  0  0  0  0\n')

            for bond_nr, bond in self.original_structure.bonds.items():
                drawn_bond = self.drawn_structure.bonds[bond_nr]
                chiral_val = None
                if (drawn_bond.atom_1.draw.is_drawn and drawn_bond.atom_2.draw.is_drawn)\
                        or (bond.atom_1.pyrrole and bond.atom_2.type == 'H')\
                        or (bond.atom_2.pyrrole and bond.atom_1.type == 'H'):

                    if drawn_bond in self.drawing.chiral_bonds:
                        wedge, atom = self.drawing.chiral_bond_to_orientation[bond]
                        if wedge == 'front':
                            chiral_val = 1
                        else:
                            chiral_val = 6
                        reverse = False
                        if atom == bond.atom_2:
                            reverse = True

                    if chiral_val:
                        if reverse:
                            molfile.write(
                                f'{str(atom_to_line_nr[bond.atom_2]).rjust(3)}{str(atom_to_line_nr[bond.atom_1]).rjust(3)}{str(self.bond_to_value[bond.type]).rjust(3)}{str(chiral_val).rjust(3)}  0  0  0\n')

                        else:
                            molfile.write(
                                f'{str(atom_to_line_nr[bond.atom_1]).rjust(3)}{str(atom_to_line_nr[bond.atom_2]).rjust(3)}{str(self.bond_to_value[bond.type]).rjust(3)}{str(chiral_val).rjust(3)}  0  0  0\n')
                    elif bond.type == 'double' and not bond.chiral and not bond.atom_1.chiral and not bond.atom_2.chiral:
                        molfile.write(
                            f'{str(atom_to_line_nr[bond.atom_1]).rjust(3)}{str(atom_to_line_nr[bond.atom_2]).rjust(3)}{str(self.bond_to_value[bond.type]).rjust(3)}  3  0  0  0\n')

                    else:
                        molfile.write(
                            f'{str(atom_to_line_nr[bond.atom_1]).rjust(3)}{str(atom_to_line_nr[bond.atom_2]).rjust(3)}{str(self.bond_to_value[bond.type]).rjust(3)}  0  0  0  0\n')

            molfile.write('M  END\n')

