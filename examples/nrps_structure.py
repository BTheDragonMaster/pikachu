from smiles2 import *
from structures import *


BONDSTRUCTURES = {
    
    'peptide': {'smiles': 'C(=O)NC',
                'structure': Smiles('C(=O)NC').smiles_to_structure(),
                'atom_1': Atom('C', 0, None, 0),
                'atom_2': Atom('N', 2, None, 0)},
    'ester_coc': {'smiles': 'C(=O)OC',
                  'structure': Smiles('C(=O)OC').smiles_to_structure(),
                  'atom_1': Atom('O', 2, None, 0),
                  'atom_2': Atom('C', 3, None, 0)}}

class NRPSStructure(Structure):
    def __init__(self, graph = None, bonds = None):
        Structure.__init__(self, graph, bonds)
        self.make_bond_lookup()
 #       self.refine_structure()

    def find_bonds(self, bond_type):
        locations = self.substructure_matching(BONDSTRUCTURES[bond_type]['structure'])
        bonds = []
        for match in locations:
            print(match)
            atom_1 = match[BONDSTRUCTURES[bond_type]['atom_1']]
            atom_2 = match[BONDSTRUCTURES[bond_type]['atom_2']]
            bond = self.bond_lookup[atom_1][atom_2]
            bonds.append(bond)

        return bonds
        
    

def main():
    string = "CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO"
 #   string = '[Co+2].[Fe+++].[O-].[OH2].[OH-].[CH4]'
 #   string = "c1ccccc1-c2ccccc2"
    string = "CCCCCCCCCC(=O)NC1C(O)C(O)C(CO)OC1Oc2c3Oc4ccc(CC5NC(=O)C(N)c6ccc(O)c(Oc7cc(O)cc(c7)C(NC5=O)C(=O)NC8C(=O)NC9C(=O)NC(C(OC%10OC(CO)C(O)C(O)C%10NC(C)=O)c%11ccc(Oc2cc8c3)c(Cl)c%11)C(=O)NC(C(O)=O)c%12cc(O)cc(OC%13OC(CO)C(O)C(O)C%13O)c%12c%14cc9ccc%14O)c6)cc4Cl"
    string = "CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
#    string = 'CC'
#    string = 'O=P([O-])([O-])[O-]'
 #   string = 'C#N'
 #   string = 'FS(F)(F)F'
    string_2 = 'c1ccccc1'
#    string = 'C(C[C@@H](C(=O))N)CN'
 #   string_2 = 'C(C[C@@H](C(=O))N)CN'
    string_2 = 'C([C@@H](C(=O)O)N)'
    string_2 = 'CC[C@H](C)[C@@H](C(=O))N'
    string_2 = 'C(=O)NC'
    smiles_1 = Smiles(string)
    smiles_2 = Smiles(string_2)

    structure_1 = smiles_1.smiles_to_structure()
    structure_1 = NRPSStructure(structure_1.graph, structure_1.bonds)
    print(structure_1.find_bonds('peptide'))
    print(structure_1.find_bonds('ester_coc'))

#    pprint(structure_1.substructure_matching(structure_2))
#    print('reverse')
#    print(structure_2.substructure_matching(structure_1))
    
 ##   structure_1.print_graph()
 #   structure_2.print_bonds()
#    for carbon in structure_1.graph:
#        print(carbon)
#        for shell in carbon.shells:
#            for orbital_set in carbon.shells[shell].orbital_sets:
#                for orbital in carbon.shells[shell].orbital_sets[orbital_set].orbitals:
#                    print(orbital)
#                    print(orbital.electrons)
 #   aromatics = structure.find_aromatic_graphs()
#    for aromatic in aromatics:
#        aromatic.print_graph()
    


if __name__ == "__main__":
    main()
