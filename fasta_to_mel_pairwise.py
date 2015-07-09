#!/usr/bin/python3

import sys
import argparse

def argument_parser():
    '''Parses the given arguments from the input'''
    parser = argparse.ArgumentParser(description='INFO ABOUT THE SCRIPT')
    parser.add_argument('input_', action='store', help='info about input')
    parser.add_argument('-l', action='store_true', dest='logs', help='write STDOUT and STDERR logs')
    arg = parser.parse_args()
    logs(arg.logs)
    return arg

def logs(logs):
    '''Checks to see if a stdout and stderr should be written, and creates logs if so'''
    if logs:
        sys.stdout = open('o.txt', 'w')
        sys.stderr = open('e.txt', 'w')

def source_file_dict(file_):
    '''takes file and turns it into a list of lists (LoL)'''
    dict_ = {}
    with open(file_, 'r') as source_file:
        for line in source_file.read().split('\n'):
            if line.startswith('>'):
                key = line 
            elif line:
                if key not in dict_.keys():
                    dict_[key] = line
                else:
                    dict_[key] += line
            else:
                break
    return dict_    

def isoform_dict(dict_):
    '''makes dict of isoforms for each species'''
    spec_dict = {}
    for key in dict_.keys():
        if key.split('; ')[-2] not in spec_dict.keys():
            spec_dict[key.split('; ')[-2]] = [[item.split('=') for item in key.split('; ')]]
        else:
            spec_dict[key.split('; ')[-2]].append([item.split('=') for item in key.split('; ')])
    gene_name = get_gene_name(spec_dict)
    return spec_dict, gene_name

def get_gene_name(dict_):
    '''gets the gene name for which the fasta contains the orthologs'''
    for list_ in dict_['species=Dmel'][0]:
        if list_[0] == 'name':
            gene_name = list_[1]
            if '-' in gene_name:
                gene_name = gene_name[:gene_name.index('-')]
            print 'Gene name =', gene_name
    return gene_name

def longest_isoform_only(dict_):
    for key, value in dict_.iteritems():
        if len(value) > 1:
            value = sorted(value, key=lambda x: float(x[-4][1]), reverse=True)[0]
            dict_[key] = '; '.join(['='.join(item) for item in value])
        else:
            dict_[key] = '; '.join(['='.join(item) for item in value[0]])
    return dict_

def remove_short_iso_fasta(fasta_dict, trim_dict):
    return {key: [value, fasta_dict[value]] for key, value in trim_dict.iteritems() if value in fasta_dict.keys()}

def write_pairwise_outputs(dict_, gene_name):
    '''Write an output file of given name'''
    keys = [key for key in dict_.keys() if key != 'species=Dmel']
    file_list = []
    for key in keys:
        output_name = '{}_{}_pairwise.fasta'.format(key[key.index('=')+2:], gene_name)
        file_list.append(output_name)
        with open(output_name, 'w+') as output_file:
            output_file.write('\n'.join(dict_['species=Dmel']))
            output_file.write('\n')
            output_file.write('\n'.join(dict_[key]))
    write_species_list(file_list)

def write_species_list(list_):
    output_name = 'species_list.txt'
    with open(output_name, 'w+') as f:
        f.write('\n'.join(list_))

def main ():
    arg = argument_parser()
    fasta_dict = source_file_dict(arg.input_)
    iso_dict, gene_name = isoform_dict(fasta_dict)
    trim_dict = longest_isoform_only(iso_dict)
    long_iso_only_dict = remove_short_iso_fasta(fasta_dict, trim_dict)
    write_pairwise_outputs(long_iso_only_dict, gene_name)

if __name__ == '__main__':
    main()
