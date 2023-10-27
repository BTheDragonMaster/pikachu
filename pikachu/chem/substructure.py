from pikachu.chem.structure import Structure

class Substructure(Structure):
    def __init__(self, graph=None, bonds=None, bond_lookup=None):
        super().__init__(graph, bonds, bond_lookup)
        self.find_cycles()