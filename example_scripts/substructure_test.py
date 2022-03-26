#!/usr/bin/env python

from pikachu.general import read_smiles, draw_smiles, highlight_subsmiles_single, highlight_subsmiles_all, highlight_subsmiles_multiple


substructure = 'OC1=CC=C(C[C@@H](N)C=O)C=C1'
substructure = "C1=CC(=CC=C1C[C@H](C(=O))N)O"
daptomycin = r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
vancomycin = r"C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O"

#read_smiles(vancomycin).print_graph()
#draw_smiles(vancomycin)

#draw_smiles(substructure)
#read_smiles(substructure).print_graph()
# highlight_subsmiles_single(substructure, daptomycin, colour='light blue', visualisation='svg', out_file='substructure_daptomycin.svg')
highlight_subsmiles_single(substructure, vancomycin, colour='light blue', visualisation='svg', out_file='substructure_vancomycin.svg')