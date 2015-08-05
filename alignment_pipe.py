#!usr/bin/python/2.7.6

import os
import sys

def list_of_proteins_and_cdss():
    '''get list of proteins and cdss'''
    return os.popen('ls *_species_prot.fasta').read().split('\n')

def run_MUSCLE_pipe(list_):
    print len(list_)
    file_name_starts = ['_'.join(line.split('_')[:2]) for line in list_ if line]
    print len(file_name_starts)
    for item in file_name_starts:
        bit = os.system('bash MUSCLE_ALIGNER.sh {0}_prot.fasta {0}_alsdjf.fasta'.format(item))
        print bit

def main():
    list_ = list_of_proteins_and_cdss()
    run_MUSCLE_pipe(list_)

if __name__ == '__main__':
    main()
