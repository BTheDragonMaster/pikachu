from pikachu import *
from reactions import BondDefiner, find_bonds, PEPTIDEBOND, ESTERCOCBOND


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


NITROGEN = GroupDefiner('carbon group', 'N', 0)
COOH = BondDefiner('carboxyl group', 'OC=O', 0, 1)

def combine_structure(structure_1, structure_2):

    for atom in structure_2.graph:
        #print(atom.nr)
        atom.nr += len(structure_1.graph)
        #print(atom.nr)

    new_graph = {}
    for atom, atoms in structure_2.graph.items():
        new_graph[atom] = atoms

    structure_2.graph = new_graph

    structure_1.make_bond_lookup()
    structure_2.make_bond_lookup()



    return structure_1, structure_2

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
    smiles_1 = Smiles("C(C(=O)O)N")
    smiles_2 = Smiles("C(C(=O)O)N")
    structure_1 = smiles_1.smiles_to_structure()
    structure_2 = smiles_2.smiles_to_structure()

    combine_structure(structure_1, structure_2)

    nitrogens = find_group(structure_1, NITROGEN)
    coohs = find_bonds(structure_2, COOH)
    nitrogen = nitrogens[0]
    cooh = coohs[0]
    structure_2.break_bond_by_nr(cooh)
    hatoms = nitrogen.neighbours
    hatom = None
    for atom in hatoms:
        if atom.type == 'H':
            hatom = atom
            break
    if hatom:
        structure_1.break_bond_between_atoms(nitrogen, hatom)

    print(structure_1.graph)
    print('\n')
    print(structure_2.graph)

