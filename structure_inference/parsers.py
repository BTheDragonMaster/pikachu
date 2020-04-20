#!/usr/bin/env python

"""
Author: Barbara Terlouw
Date: 05/07/2018

Parsing functions for general use

"""

import json
import pprint

def parse_mass_file(mass_file):
    mass_file.readline()

    for line in mass_file:
        line = line.strip()
        if line:
            compound, mono_iso, formula, organism  = line.split('\t')
            compound = compound.strip()
            formula = formula.strip()
            organism = organism.strip()
            mono_iso = float(mono_iso.replace(',', '.'))
            yield compound, mono_iso, formula, organism

def parse_info_file(info_file):
    info_file.readline()

    for line in info_file:
        acc, nr, contig, in_ann, start, end, prot_ID, gene_ID, gene_func,\
             cat, tailor, evidence, KO, publ, comm, target_aa, strand = \
             line.split('\t')
        nr = int(nr)
        start = int(start)
        end = int(end)

        yield acc, nr, contig, in_ann, start, end, prot_ID, gene_ID, \
              gene_func, cat, tailor, evidence, KO, publ, comm, target_aa, \
              strand
        
def parse_compound_file(compound_file):
    
    for line in compound_file:
        line = line.strip()
        if line:
            split_line = line.split('\t')
            BGC = split_line[0]
            compounds = split_line[1:]
            yield BGC, compounds
        
    
def parse_sandpuma_output(sandpuma_file):
    for line in sandpuma_file:
        line = line.strip()
        if line:
            ID, pred_mono, pred_snn, snn_score, asm, svm, phmm, pid,\
                ensemble, sandpuma = line.split('\t')
            yield ID, pred_mono, pred_snn, snn_score, asm, svm, phmm, pid,\
                  ensemble, sandpuma

def parse_smiles_file(smiles_file):
    for line in smiles_file:
        line = line.strip()
        ID, smiles = line.split('\t')
        yield ID, smiles

def parse_subsmiles_file(smiles_file):
    for line in smiles_file:
        line = line.strip()
        if line:
            ID = line.split('\t')[0]
            smiles = line.split('\t')[1:]
            yield ID, smiles

def parse_mol_file(mol_file):
    atom_dict = {}
    atom_nr = 0
    structure_graph = {}

    for line in mol_file:
        line_stripped = line.strip()
        
        if line:
            if line_stripped.startswith(tuple('0123456789-')):
                data = line_stripped.split()
                if len(data) == 11:
                    pass

                elif len(data) == 16:
                    atom_nr += 1
                    atom_dict[atom_nr] = data[3]

                elif len(data) <= 4:
                    nr_1 = line[:3]
                    nr_2 = line[3:6]
                    nr_3 = line[6:9]
                    nr_4 = line[9:12]
                    
                    atom_nr_1 = int(nr_1.strip())
                    atom_nr_2 = int(nr_2.strip())
                    atom_1 = (atom_dict[atom_nr_1], atom_nr_1)
                    atom_2 = (atom_dict[atom_nr_2], atom_nr_2)
                    bond_nr = int(nr_3.strip())
                    
                    if atom_1 in structure_graph:
                        structure_graph[atom_1] += bond_nr * [atom_2]
                    else:
                        structure_graph[atom_1] = bond_nr * [atom_2]

                    if atom_2 in structure_graph:
                        structure_graph[atom_2] += bond_nr * [atom_1]
                    else:
                        structure_graph[atom_2] = bond_nr * [atom_1]
                else:
                    pass

    return structure_graph
                
                
        
def parse_monomers(monomers_file):
    monomers_file.readline()
    for line in monomers_file:
        line = line.strip()
        if line:
            record = line.split('\t')
            yield record
            

def make_clan_dict(clan_file):
    """Return dict of {PFAM clan: [PFAM family, ->], ->} from txt file

    Input:
    clan_file: file, two columns with the first containing clan accesions, and
        the second containing family accessions

    Output:
    clan_dict: dict of {PFAM clan: [PFAM family, ->], ->}, with PFAM clan and
        pfam family both str
    """
    clan_dict = {}
    
    for line in clan_file:
        line = line.strip()
        clan, family = line.split()
        if clan in clan_dict:
            
            if family not in clan_dict[clan]:
                clan_dict[clan] += [family]
        else:
            clan_dict[clan] = [family]

    return clan_dict

def add_entry_structure_dict(json_data, BGC, structure_dict):

    structure_dict[BGC] = []
    for compound in json_data["compounds"]:
        structure_dict[BGC] += [compound["compound"].lower()]
    

def parse_json(json_file):
    """Return json dict from json file extracted from MIBiG database

    Input:
    json_file: file, .json file from MIBiG database

    Output:
    json_data: dict, quite extensive..
    """
    json_data = json.load(json_file)["general_params"]
    return json_data

def add_entry_BGC_dict(json_data, BGC_dict, BGC_name):
    """Add entry to BGC dict

    Input:
    json_data: dict, quite extensive..
    BGC_dict: dict of {BGC_name: [biosynthetic class, ->], ->}
    BGC_name: str, name of BGC cluster in the form BGC%07d % BGC number
    """
    BGC_dict[BGC_name] = get_biosyn_classes(json_data)
    

def get_biosyn_classes(json_data):
    """Return biosynthetic gene cluster classes from json data

    Input:
    json_data: dict, quite extensive..

    Output:
    biosyn_classes: list of str, with each str a biosynthetic gene clusrer
        class
    """
    biosyn_classes = json_data["biosyn_class"]
    return biosyn_classes

def get_smiles(json_data):
    """Return smiles strings from json data

    Input:
    json_data: dict, quite extensive..

    Output:
    smiles: list of str, with each str a smiles of a structure produced by the
        biosynthetic gene cluster described in the .json file
    """

    smiles = []
    compounds = []

    for structure in json_data["compounds"]:
        try:
            smile = structure["chem_struct"]
            compound = structure["compound"]
            if smile:
                smiles += [smile]
                compounds += [compound]
        except KeyError:
            pass

    return smiles, compounds
                    

def make_reverse_dict(some_dict):
    """Return reverse dict from dict, such that entries point to keys

    Input:
    some_dict: dict of {key: [new_key, ->], ->}

    Output:
    reverse_dict: dict of {new_key: [key, ->, ->}
    """
    reverse_dict = {}
    for key in some_dict:
        new_keys = some_dict[key]
        for new_key in new_keys:
            if new_key in reverse_dict:
                
                reverse_dict[new_key] += [key]
            else:
                reverse_dict[new_key] = [key]
    return reverse_dict
        

def parse_fasta(fasta_file):
    """Yield gene/protein IDs and sequences from a .fasta file

    Input:
    fasta_file: file, .fasta file

    Output:
    record: list of [ID, sequence], with ID and sequence both str.
    """
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        
        if line.startswith(">"):
            if sequence:
                record = [ID, sequence]
                yield record
                
                ID = line[1:]
                sequence = ''
            else:
                ID = line[1:]
                
        else:
            sequence += line

    record = [ID, sequence]
    yield record

def parse_protein_ID(fasta_ID):
    """Return protein ID from MIBiG fasta ID

    Input:
    fasta_ID: str, ID taken directly from MIBiG .fasta file

    Output:
    protein_ID: str, protein ID
    """
    fasta_ID = fasta_ID.strip()
    protein_ID = fasta_ID.split("|")[-1].strip()
    return protein_ID

def parse_ordered_pfam(ordered_file):
    first_line = True
    record = []
    
    for line in ordered_file:
        line = line.strip()
        if line.startswith("Clan:"):
            
            if first_line:
                first_line = False
                ID = line.split(": ")[1]
                record += [ID]
            else:
                yield record
                record = []
                ID = line.split(": ")[1]
                record += [ID]
            
        else:
            record += [tuple(line.split('\t'))]
    yield record

def parse_pfam_sorted(ordered_pfam):
    for line in ordered_pfam:
        line = line.strip()
        if line.startswith("PF"):
            record = line.split("\t")
            yield record
                
def add_entry_to_clan_pfam_dict(record, clan_pfam_dict):
    ID, prot_info = record
    if ID in clan_pfam_dict:
        clan_pfam_dict[ID] += [prot_info]
    else:
        clan_pfam_dict[ID] = [prot_info]
    
            

def parse_cluster_ID(fasta_ID):
    """Return BGC ID from MIBiG fasta ID

    Input:
    fasta_ID: str, ID taken directly from MIBiG .fasta file

    Output:
    cluster_ID: str, BGC ID
    """
    fasta_ID = fasta_ID.strip()
    cluster_ID = fasta_ID.split("|")[0].strip()
    return cluster_ID


def parse_pfam_domtbl(domtbl_file):
    
    for line in domtbl_file:
        if line.startswith("#"):
            continue
        else:
            line = line.strip()
            target, target_acc, target_l, query, query_acc, query_l, f_eval,\
                    f_score, f_bias, current_dom, total_dom, c_eval, i_eval,\
                    d_score, d_bias, hmm_start, hmm_stop, seq_start, seq_stop,\
                    env_start, env_stop, reliability = line.split()[:22]
            target_description = line.split()[22:]
            target_description = ' '.join(target_description)
            
            protein = parse_protein_ID(query)
            cluster = parse_cluster_ID(query)

            target_l = int(target_l)
            query_l = int(query_l)

            current_dom = int(current_dom)
            total_dom = int(total_dom)
            
            f_eval = float(f_eval)
            c_eval = float(c_eval)
            i_eval = float(i_eval)

            f_score = float(f_score)
            d_score = float(d_score)

            f_bias = float(f_bias)
            d_bias = float(d_bias)

            hmm_start = int(hmm_start)
            hmm_stop = int(hmm_stop)

            seq_start = int(seq_start)
            seq_stop = int(seq_stop)

            env_start = int(env_start)
            env_stop = int(env_stop)

            reliability = float(reliability)

            yield target, target_acc, target_l, query, query_acc, query_l, f_eval,\
                  f_score, f_bias, current_dom, total_dom, c_eval, i_eval,\
                  d_score, d_bias, hmm_start, hmm_stop, seq_start, seq_stop,\
                  env_start, env_stop, reliability, target_description,\
                  cluster, protein
            

            
            

            

            
            
                
            
    
