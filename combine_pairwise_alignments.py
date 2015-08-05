#!/usr/bin/python3

import sys
import argparse
import copy

def argument_parser():
    '''Parses the given arguments from the input'''
    parser = argparse.ArgumentParser(description='Combine pairwise alignments of mel gene against other species orthologs')
    parser.add_argument('input_', nargs='+', action='store', help='Pairwise alignment')
    parser.add_argument('-ll', action='store', dest='line_length', type=int, default=60, help='Indicate length of lines in alignment file (default = 60)')
    parser.add_argument('-nuc', action='store', dest='nuc', help='substitute nucleotides in for prot')
    parser.add_argument('-l', action='store_true', dest='logs', help='write STDOUT log')
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
        dict_, gene_name = make_dict(file_, dict_)
    return dict_, gene_name

def make_dict(file_, dict_):
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
    print '\nThese are the insertion indicies for each species'
    old_counter = 0
    for key, value in mel_seq.iteritems():
        if 'mel_final' not in locals():
            mel_final = [char for char in value if char != '-']

        temp_value, temp_mel_final, new_mel_final = value[:], mel_final[:], []
        while len(temp_value) > 0 and len(temp_mel_final) > 0:
            for i, (a, b) in enumerate(zip(temp_value, temp_mel_final)):
                if a == b:
                    new_mel_final.append(a)
                    temp_value.pop(i)
                    temp_mel_final.pop(i)
                    break
                elif a == '-':
                    new_mel_final.append(a)
                    temp_value.pop(i)
                    break
                elif b == '-':
                    new_mel_final.append(b)
                    temp_mel_final.pop(i)
                    break
                else:
                    print 'Whoops, these should match and done', a, '---', b
                    raw_input()
        mel_final = new_mel_final
    return mel_final

def mk_insertion_dict(mel_seq_original, mel_final):
    gaps = [i for i, gap in enumerate(mel_final) if gap == '-']
    insert_dict = {}
    for key, value in mel_seq_original.iteritems():
        insert_dict[key] = []
        for i, nucl in enumerate(value):
            if i in gaps and nucl != '-':
                value.insert(i, '-')
            elif nucl == '-':
                insert_dict[key].append(i) 
        print key, ', '.join(map(str, insert_dict[key]))
    print 'Dmel length', len(mel_final)
    return ''.join(mel_final), insert_dict

def add_gaps_to_species(spec_seq, mel_w_gaps, insert_dict):
    '''Add gaps to each of the species genes'''
    gaps = [i for i, gap in enumerate(mel_w_gaps) if gap == '-']
    for key, value in spec_seq.iteritems():
        for index in gaps:
            if index not in insert_dict[key]:
                value.insert(index, '-') 
    spec_seq['mel'] = list(mel_w_gaps)
    return spec_seq

def substitute_nucleotides(spec_prot_align, spec_nucl_seq):
    if spec_nucl_seq:
        nuc_dict = {}
        try:
            with open(spec_nucl_seq, 'r') as source_file:
                gene_name = source_file.name.split('_')[1]
                for i, line in enumerate(source_file.read().split('\n')):
                    if i == 0:
                        key, value = line[line.index('species=D') + 9:line.index('species=D') + 12], ''
                    elif line.startswith('>'):
                        if len(spec_prot_align[key]) - spec_prot_align[key].count('-') + 1 == len(value) / 3:
                            n = 3
                            value = [value[i:i+n] for i in range(0, len(value), n)]
                            nuc_dict[key] = value
                        key, value = line[line.index('species=D') + 9:line.index('species=D') + 12], ''
                    elif line:
                        value += line
                else:
                    if len(spec_prot_align[key]) - spec_prot_align[key].count('-') + 1 == len(value) / 3:
                        n = 3
                        value = [value[i:i+n] for i in range(0, len(value), n)]
                        nuc_dict[key] = value

            for key, value in spec_prot_align.iteritems():
                nucl_seq = nuc_dict[key]
                new_value, count = '', 0
                for i, amino_acid in enumerate(value):
                    if amino_acid == '-':
                        new_value += '---'
                    elif amino_acid != '-':
                        new_value += nucl_seq[count]
                        count += 1
                spec_prot_align[key] = list(new_value)
            return spec_prot_align
        except IOError:
            return spec_prot_align
    else:
        return spec_prot_align

def fill_out_values(gapped_seq):
    max_len = max([len(value) for value in gapped_seq.values()])
    for key, value in gapped_seq.iteritems():
        gapped_seq[key] = value + ['-'] * (max_len - len(value))
    return gapped_seq

def check_per_align(gapped_seq, mel):
    '''Checks to see if there is perfect alignment of the amino acid'''
    perf_align = []
    values = gapped_seq.values()
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
    keys = sorted([key for key in spec_seq.keys() if key != 'mel'])
    keys.insert(0, 'mel')
    keys.append('perf')
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
    mel_w_gaps = combine_mel_alignments(mel_seq)
    mel_w_gaps, insert_dict = mk_insertion_dict(mel_seq, mel_w_gaps)
    gapped_seq = add_gaps_to_species(spec_seq, mel_w_gaps, insert_dict)
    gapped_seq = substitute_nucleotides(gapped_seq, arg.nuc)
    gapped_seq = fill_out_values(gapped_seq)
    perf_align = check_per_align(gapped_seq, mel_w_gaps)
    write_output(mel_w_gaps, gapped_seq, arg.line_length, perf_align, gene_name)

if __name__ == '__main__':
    main()
