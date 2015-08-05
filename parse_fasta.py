#!/usr/bin/python3

import sys
import argparse

def logs(logs):
    '''Checks to see if a stdout and stderr should be written, and creates logs if so'''
    if logs:
        sys.stdout = open('{}_log.txt'.format(sys.argv[0][:-3]), 'w')

def source_file_input(input_):
    '''
    takes fasta file and creates dict with FBgn being key and sequence being value
    '''
    with open(input_, 'r') as source_file:
        dict_of_set_species = {}
        list_ = [line for line in source_file.read().split('>') if line]
        if '#' in list_[-1]:
            list_[-1] = list_[-1][:list_[-1].index('#')]

        for char in list_[0][-10:].strip():
            if char not in ['A', 'T', 'G', 'C', 'N']:
                type_ = 'prot'
                break
        else:
            type_ = 'nucl'

        for line in list_:
            specie = line[line.index('species=') + 8 : line.index('species=') + 12]
            length = line[line.index('length=') + 7 :]
            length = length[:length.index(';')]
            if not line.startswith('FBgn'):
                key = line[line.index('parent=') + 7 : line.index('parent=') + 18]
            elif line.startswith('FBgn'):
                key = line[:12]
            try:
                dict_of_set_species[key].append((specie,length, line))
            except KeyError:
                dict_of_set_species[key] = [(specie, length, line)]
    return dict_of_set_species, type_

def orthologs_file_input(orthos):
    '''opens ortholog file and makes it into a dict'''
    with open(orthos, 'r') as source_file:
        return {line.split('\t')[0]: line.split('\t')[1:] for line in source_file.read().split('\n')[1:]}

def make_gene_dict(species_dict, ortho_dict):
    gene_dict = {}
    for key, value in ortho_dict.iteritems():
        for ortho_FBgn in value:
            if ortho_FBgn in species_dict.keys():
                try:
                    gene_dict[key] += species_dict[ortho_FBgn]
                except KeyError:
                    gene_dict[key] = species_dict[ortho_FBgn]
    return gene_dict

def orthologs_present(dict_):
    '''checks to see if the dictionary value contains more than just D.mel variants'''
    ortho_dict = {}
    for key, value in dict_.iteritems():
        species = list(set([bit[0] for bit in value]))
        lengths = [bit[:2] for bit in value]
        if len(species) != 1:
            ortho_dict[key] = sorted(value, key = lambda x : x[1], reverse=True)
    return ortho_dict

def longest_iso_only(dict_):
    final_dict = {}
    for key, value in dict_.iteritems():
        species_included = []
        longest_only = []
        for item in value:
            if item[0] not in species_included:
                species_included.append(item[0])
                longest_only.append(item[2])
        final_dict[key] = longest_only
    return final_dict

def write_output_file(final_dict, type_):
    '''Write an output file of given name'''
    for key, value in final_dict.iteritems():
        output_name = '{}_species_{}.fasta'.format(key, type_)
        with open(output_name, 'w') as output_file:
            value.insert(0, '')
            output = '>'.join(value)
            output_file.write(output)

def main(arg):
    logs(arg.logs)
    dict_of_set_species, type_ = source_file_input(arg.input_)
    ortho_dict = orthologs_file_input(arg.orthos)
    gene_dict = make_gene_dict(dict_of_set_species, ortho_dict)
    ortho_dict = orthologs_present(gene_dict)
    final_dict = longest_iso_only(ortho_dict)
    write_output_file(final_dict, type_)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parses fasta containing many different genes to many fastas, each containing the longest ortholog variant for the gene')
    parser.add_argument(
            'input_', 
            action='store', 
            help='info about input'
            )
    parser.add_argument(
            'orthos', 
            action='store', 
            help='ortholog_file'
            )
    parser.add_argument(
            '-l', 
            action='store_true', 
            dest='logs', 
            help='write STDOUT log'
            )
    arg = parser.parse_args()
    main(arg)
