import sys
import argparse

def make_ortho_dict(ortho_input):
    '''makes dict out of the final ortholog list'''
    ortho_dict = {}
    with open(ortho_input, 'r') as f:
        for line in f.read().split('\n'):
            if line and not line.startswith('#'):
                try:
                    ortho_dict[line.split('\t')[0]].append(line)
                except KeyError:
                    ortho_dict[line.split('\t')[0]] = [line]
    return ortho_dict

def open_FBgn_query_list(FBgn_input):
    '''
    Opens file containing list of FBgn's that you want to match against the ortholog list,
    and converts into list.  For Flybase searches, use quicksearch (eg. olfactory receptors), 
    then on genes list, click HitList conversion tools and pick FlyBase IDs file option
    '''
    with open(FBgn_input, 'r') as f:
        return [line for line in f.read().split('\n') if line.strip()]

def build_query(FBgn_list, ortho_dict):
    '''Makes list containing D.mel FBgn and the orthologs from the species defined in the species_list'''
    species_list = ['Dere', 'Dsec', 'Dvir']
    for FBgn in FBgn_list[:]:
        try:
            for line in ortho_dict[FBgn]:
                for specie in species_list:
                    if specie in line:
                        FBgn_list.append(line.split('\t')[5])
        except KeyError:
            print FBgn
    return FBgn_list

def write_FBgn_output_file(FBgn_list, output_name):
    '''Writes D.mel and orthologs FBgn to output file'''
    if output_name == None:
        output_name = 'query_list'
    with open('{}.txt'.format(output_name), 'w') as f:
        f.write('\n'.join(FBgn_list))

def main(arg):
    ortho_dict = make_ortho_dict(arg.ortho)
    FBgn_list = open_FBgn_query_list(arg.input_)
    FBgn_list = build_query(FBgn_list, ortho_dict)
    write_FBgn_output_file(FBgn_list, arg.out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Makes list of D.mel FBgn and orthologs from a list of only D.mel FBgns')
    parser.add_argument(
        'input_', 
        action='store', 
        help='D.mel FBgn list of genes you want orthologs of')
    parser.add_argument(
        'ortho', 
        action='store', 
        help='Table containing gene name, D.mel FBgn, and identified orthologs')
    parser.add_argument(
        '-o', 
        action='store', 
        dest='out', 
        help='Define output name, otherwise will be query_list.txt')
    parser.add_argument(
        '-l', 
        action='store_true', 
        dest='logs', 
        help='write STDOUT')
    arg = parser.parse_args()

    main(arg)
