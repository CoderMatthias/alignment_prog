#!/usr/bin/python3

import sys
import argparse

def argument_parser():
    '''Parses the given arguments from the input'''
    parser = argparse.ArgumentParser(description='INFO ABOUT THE SCRIPT')
    parser.add_argument('input_', action='store', help='info about input')
    parser.add_argument('-l', action='store_true', dest='logs', help='write STDOUT and STDERR logs')
    arg = parser.parse_args()
    logs(spec_name, arg.logs)
    return arg, spec_name

def logs(name, logs):
    '''Checks to see if a stdout and stderr should be written, and creates logs if so'''
    if logs:
        sys.stdout = open('8_{}_synteny_o.tsv'.format(name), 'w')
        sys.stderr = open('8_{}_synteny_e.tsv'.format(name), 'w')

def source_file_LoL (sys_argv):
    '''takes file and turns it into a list of lists (LoL)'''
    with open(sys_argv, 'r') as source_file:
        header = source_file.readline().strip().split('\t')
        line_list = [line for line in  source_file.read().split('\n') if line.strip() and not line.startswith('#')]
    return [line.split('\t') for line in line_list], header

def write_synteny_output (synteny_file, header, species_name):
    '''Write an output file of given name'''
    output_name = '8_{}_synteny_analyzed.tsv'.format(species_name) 
    print 'Output saved as:', output_name
    with open(output_name, 'w+') as output_file:
        output_file.write('{}\n'.format('\t'.join(header)))
        for line in synteny_file:
            output_line = '{}\n'.format('\t'.join(line))
            output_file.write(output_line)

def main ():
    arg, species_name = argument_parser()

if __name__ == '__main__':
    main()
