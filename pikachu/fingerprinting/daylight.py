import hashlib

from pikachu.chem.atom_properties import ATOM_PROPERTIES


class Daylight:
    def __init__(self, atom, structure):
        self.atom = atom
        self.structure = structure
        self.d1 = 0
        self.d2 = 0
        self.d3 = 0
        self.d4 = 0
        self.d5 = 0
        self.d6 = 0
        self.d7 = 0
        self.daylight = tuple()
        self.set_daylight_properties()

    def set_daylight_properties(self):
        self.d1 = self.get_heavy_neighbours()
        self.d2 = self.get_valence_minus_h()
        self.d3 = ATOM_PROPERTIES.element_to_atomic_nr[self.atom.type]
        self.d4 = ATOM_PROPERTIES.element_to_amu[self.atom.type]
        self.d5 = self.atom.charge
        self.d6 = self.get_hydrogen_number()
        self.d7 = self.atom_in_cycle()
        self.daylight = (self.d1, self.d2, self.d3, self.d4, self.d5, self.d6, self.d7)

    def get_hash(self):
        daylight_hash = hashlib.sha256()
        for attribute in self.daylight:
            daylight_hash.update(str(attribute).encode())

        # return int.from_bytes(hashlib.sha256(b"H").digest()[:4], 'little')
        #print(hash(tuple(self.daylight)))
       # print(hash(tuple(self.daylight)))
        return hash(tuple(self.daylight))

    def atom_in_cycle(self):
        cycles = self.structure.cycles.all_cycles

        for cycle in cycles:
            if self.atom in cycle:
                return 1

        return 0

    def get_heavy_neighbours(self):
        heavy_neighbours = 0
        for atom in self.atom.neighbours:
            if atom.type != 'H':
                heavy_neighbours += 1
        return heavy_neighbours

    def get_hydrogen_number(self):
        hydrogen_nr = 0
        for atom in self.atom.neighbours:
            if atom.type == 'H':
                hydrogen_nr += 1
        return hydrogen_nr

    def get_valence_minus_h(self):
        valence = self.atom.calc_bond_nr()

        return valence - self.get_hydrogen_number()