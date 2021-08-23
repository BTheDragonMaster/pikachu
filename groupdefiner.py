from pikachu import *

class GroupDefiner():
    def __init__(self, name, smiles, atom_nr_1):
        self.name = name
        self.smiles = Smiles(smiles)
        self.structure = self.smiles.smiles_to_structure()
        self.find_atom(atom_nr_1)

    def __repr__(self):
        return self.name

    def find_atom(self, atom_nr_1):

        self.atom_1 = None

        for atom in self.structure.graph:
            if atom.nr == atom_nr_1:
                self.atom_1 = atom

            if self.atom_1:
                break

        if not self.atom_1:
            raise Exception("Can't find atoms adjacent to bond.")



CARBON = GroupDefiner('carbon group', 'C', 0)
NITROGEN = GroupDefiner('nitrogen group', 'N', 0)
CARBOXYLIC_ACID = GroupDefiner('carboxylic acid group', 'OC=O', 1)
AMP_PHOSPHOR = GroupDefiner('AMP_PHOSPHOR', 'COP(=O)(O)', 2)
PPi_O = GroupDefiner('PPi_O', 'OP(=O)(O)OP(=O)(O)O', 0)
HYDROXYL = GroupDefiner('HYDROXYL', 'O', 0)
ESTER = GroupDefiner('ESTER', 'C(=O)O', 0)

def find_group(structure, group_type):
    locations = structure.find_substructures(group_type.structure)
    groups = []
    for match in locations:
        group = match[group_type.atom_1]
        groups.append(group)

    return groups



if __name__ == "__main__":
    string = "CCCCCCCCCC(=O)NC1C(O)C(O)C(CO)OC1Oc2c3Oc4ccc(CC5NC(=O)C(N)c6ccc(O)c(Oc7cc(O)cc(c7)C(NC5=O)C(=O)NC8C(=O)NC9C(=O)NC(C(OC%10OC(CO)C(O)C(O)C%10NC(C)=O)c%11ccc(Oc2cc8c3)c(Cl)c%11)C(=O)NC(C(O)=O)c%12cc(O)cc(OC%13OC(CO)C(O)C(O)C%13O)c%12c%14cc9ccc%14O)c6)cc4Cl"
    string = "CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
    smiles_1 = Smiles(string)
    structure_1 = smiles_1.smiles_to_structure()

    carbon_group = find_group(structure_1, CARBON)
    print(carbon_group)
    carboxylic = find_group(structure_1, CARBOXYLIC_ACID)
    print(carboxylic)
