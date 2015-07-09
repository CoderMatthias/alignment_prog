#!/usr/bin/python3

import sys
import argparse
import time

def argument_parser():
    '''Parses the given arguments from the input'''
    parser = argparse.ArgumentParser(description='Combine pairwise alignments of mel gene against other species orthologs')
    parser.add_argument('input_', nargs='+', action='store', help='Pairwise alignment')
    parser.add_argument('-a', action='store', dest='align_file', help='File containing alignments')
    parser.add_argument('-ll', action='store', dest='line_length', type=int, default=60, help='Indicate length of lines in alignment file (default = 60)')
    parser.add_argument('-l', action='store_true', dest='logs', help='write STDOUT and STDERR logs')
    arg = parser.parse_args()
    logs(arg.logs)
    return arg

def logs(logs):
    '''Checks to see if a stdout and stderr should be written, and creates logs if so'''
    if logs:
        sys.stdout = open('{}_o.tsv'.format(sys.argv[0]), 'w')

def source_file_dict(input_):
    '''
    makes dict where key=species name, 
    value[0][0]=alignment info line for mel,  value[0][1]=mel aligned sequence, 
    value[1][0]=alignment info line for spec, value[1][1]=spec aligned sequence
    '''
    dict_ = {}
    for file_ in input_:
        with open(file_, 'r') as source_file:
            species_name = source_file.name.split('_')[0]
            gene_name = source_file.name.split('_')[1]
            for i, line in enumerate(source_file.read().split('\n')):
                if i == 0:
                    value = [line]
                elif line.startswith('>'):
                    dict_[species_name] = [value]
                    value = [line]
                elif line:
                    if len(value) == 1:
                        value.append(line)
                    else:
                        value[1] += line
            else:
                dict_[species_name].append(value)
    return dict_, gene_name

def make_mel_and_species_dicts(dict_):
    '''make two dicts: one containing mel portion of alignment, other containing spec portion of alignment'''
    mel_seq, spec_seq, insert_dict = {}, {}, {}
    for key, value in dict_.iteritems():
        mel_seq[key] = list(value[0][1])
        spec_seq[key] = list(value[1][1])
    return mel_seq, spec_seq
   
def combine_mel_alignments(mel_seq):
    '''Make a mel gene containing all gaps from each pairwise alignment'''
    insert_dict = {}
    temp = []
    dict1 = mel_seq.copy()
    print '\nThese are the insertion indicies for each species'
    for key, value in mel_seq.iteritems():
        insert_dict[key] = []
        if 'mel_final' not in locals():
            mel_final = [char for char in value if char != '-']
### THIS IS FOR ADDING - TO THE SPECIES MEL PORTIONS AND NOTING INSERTIONS LOCATIONS ###
        gaps = [i for i, gap in enumerate(mel_final) if gap == '-']
        for index in gaps:
            for i, nucl in enumerate(value):
                if i == index and nucl != '-':
                    value.insert(i, '-')
                elif i == index and nucl == '-':
                    insert_dict[key].append(index)
        i_mel_final = []
        count = 0
        for i, nucl in enumerate(value):
            if nucl == '-' and i not in gaps:
                i_mel_final.append(i)
                for value in insert_dict.values():
                    for ind, ite in enumerate(value):
                        if i < ite:
                            value[ind] = int(ite) + 1
            elif nucl == '-' and i in gaps:
                count += 1
        insert_dict[key] = insert_dict[key] + i_mel_final
        print key, ', '.join(map(str, insert_dict[key]))
        for i, index in enumerate(i_mel_final):
            mel_final.insert(index, '-')
    gaps = [i for i, gap in enumerate(mel_final) if gap == '-']
    for key, value in dict1.iteritems():
        print 'heyo', key
        tempt = []
        for i, nucl in enumerate(value):
            if i in gaps and nucl != '-':
                print 'no gap', i
                value.insert(i, '-')
            elif i in gaps and nucl == '-':
                print 'gap', i
                tempt.append(i)
            elif nucl == '-':
                print 'gap', i
                tempt.append(i)
        print ', '.join(map(str, tempt))

    print 'Dmel length', len(mel_final)
    return ''.join(mel_final), insert_dict

def add_gaps_to_species(spec_seq, mel_w_gaps, insert_dict):
    '''Add gaps to each of the species genes'''
    gaps = [i for i, gap in enumerate(mel_w_gaps) if gap == '-']
    for key, value in spec_seq.iteritems():
        for index in gaps:
            if index not in insert_dict[key]:
                value.insert(index, '-') 
    return spec_seq

def check_per_align(gapped_seq, mel):
    perf_align = []
    values = gapped_seq.values() + [mel]
    for i in range(len(values[0])):
        temp = []
        for spec in values:
            temp.append(spec[i])
        if len(set(temp)) == 1:
            perf_align.append('*')
        else:
            perf_align.append(' ')
    return ''.join(perf_align)

def write_output(mel_w_gaps, spec_seq, per_line, perf_align, gene_name):
    '''Write an output file'''
    keys = sorted(spec_seq)
    keys.insert(0, 'mel')
    keys.append('perf')
    spec_seq['mel'] = mel_w_gaps
    spec_seq['perf'] = perf_align
    output_name = '{}_alignment.txt'.format(gene_name)
    print '\nOutput saved as:', output_name
    print
    with open(output_name, 'w+') as output_file:
        writing = True
        while writing == True:
            for key in keys:
                output_line = '{}\t'.format(key)
                if spec_seq[key] == 'break':
                    writing = False
                    break
                elif len(spec_seq[key]) > per_line:
                    output_line += ''.join(spec_seq[key])[:per_line]
                    spec_seq[key] = ''.join(spec_seq[key])[per_line:]
                elif len(spec_seq[key]) <= per_line:
                    output_line += spec_seq[key]
                    spec_seq[key] = 'break'
                output_file.write('{}\n'.format(output_line))
            output_file.write('\n')

def main ():
    arg = argument_parser()
    full_dict, gene_name = source_file_dict(arg.input_)
    mel_seq, spec_seq = make_mel_and_species_dicts(full_dict)
    mel_w_gaps, insert_dict = combine_mel_alignments(mel_seq)
    gapped_seq = add_gaps_to_species(spec_seq, mel_w_gaps, insert_dict)
    perf_align = check_per_align(gapped_seq, mel_w_gaps)
    write_output(mel_w_gaps, gapped_seq, arg.line_length, perf_align, gene_name)

if __name__ == '__main__':
    main()
