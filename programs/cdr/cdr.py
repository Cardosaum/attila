'''
cdr.py - Analyze an "aafreq.txt" input file and
         return all features of CDR3 sequences

Usage:
       python3 cdr.py <INPUT_FILE> <OUTPUT_DIRECTORY>

Copyright:  (c)  2019  Matheus Cardoso  <https://github.com/cardosaum>

License:  Apache 2.0  <https://www.apache.org/licenses/LICENSE-2.0>
'''


import re
import sys
import timeit
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import logging
import csv
import collections
import pathlib


logger = logging.getLogger('cdrlog')
hdlr = logging.FileHandler('cdrlog.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)


def extract_cdr3(file):
    '''
    Read all input file and retur a dictionary in format:
            'some_cdr3': <quantity_x>
            'other_cdr3': <quantity_y>
    Where <quantity> stands for the number of repeated <cdr3> pattern
    '''

    cdr3_regex = re.compile(r'((.+)(C)(.+)(C)(.{2})(.+)(WG.G)(.+)?)')
    dic_cdr3 = {}
    # TODO: Attila is generating files with non UTF-8 characters, like 'ÿ'
    with open(file, encoding="ISO-8859-1") as f:
        read = False
        for line in f:
            if line.startswith('>'):

                # That way, the next line we know that is a cdr sequence
                read = True

            elif read:
                parse = cdr3_regex.search(line)

                # Correspond to `(.+)` in midle of `(.{2})` and `(WG.G)`
                # in `cdr3_regex` variable
                cdr3 = parse[7]

                if cdr3 not in dic_cdr3:
                    dic_cdr3.setdefault(cdr3, 1)

                else:
                    dic_cdr3[cdr3] += 1

                # Need to be avaluated in next loop
                read = False
    dic_cdr3 = sorted(dic_cdr3.items(), key=lambda kv: kv[1], reverse=True)
    dic_cdr3 = collections.OrderedDict(dic_cdr3)
    return dic_cdr3


def aa_groups(sequence):

    group = {
            'aliphatic': ['G', 'A', 'P', 'V', 'L', 'I', 'M'],
            'aromatic': ['F', 'Y', 'W'],
            'neutral': ['S', 'T', 'C', 'N', 'Q'],
            'positive': ['K', 'H', 'R'],
            'negative': ['D', 'E'],
            }

    result = {
            'aliphatic': 0,
            'aromatic': 0,
            'neutral': 0,
            'positive': 0,
            'negative': 0,
            'invalid': 0,
            }

    not_listed = set()

    for aa in sequence.upper():
        if aa in group['aliphatic']:
            result['aliphatic'] += 1
        elif aa in group['aromatic']:
            result['aromatic'] += 1
        elif aa in group['neutral']:
            result['neutral'] += 1
        elif aa in group['positive']:
            result['positive'] += 1
        elif aa in group['negative']:
            result['negative'] += 1
        else:
            not_listed.add(aa)
            result['invalid'] += 1

    return result


# Write to a file the attributes of CDR3 sequence
# In order, the attributes written are:
# TODO: Place HEADER here

def write_cdr3_attributes(file, output, prefix='', suffix=''):
    '''
Take <file> as absolute path to the file that will be analised
and <output> as absolute path do the directory where the result analisys
will be save.

Hence,

<file>: absolute path to analysed file
<output>: absolute path to output file
    '''

    # Write HEADER to file
    output_file = pathlib.Path(output).resolve()

    if not output_file.parent.is_dir():
        raise EOFError('<output> variable must be a valid directory')

    header = ['cdr3', 'quantity', 'length', 'MW', 'AV', 'IP', 'flex', 'gravy',
              'SSF_Helix', 'SSF_Turn', 'SSF_Sheet', 'n_A', 'n_C', 'n_D', 'n_E',
              'n_F', 'n_G', 'n_H', 'n_I', 'n_K', 'n_L', 'n_M', 'n_N', 'n_P',
              'n_Q', 'n_R', 'n_S', 'n_T', 'n_V', 'n_W', 'n_Y', 'aliphatic',
              'aromatic', 'neutral', 'positive', 'negative', 'invalid', 'file']

    aa_error = 0
    aa_total = 0

    with output_file.open("w") as out:
        out = csv.writer(out)
        out.writerow(header)
        for cdr3, quantity in extract_cdr3(file).items():
            try:
                attributes = []
                # print(f'CDR3:\t{cdr3}'.ljust(60)+f'Quantity: {quantity}'.rjust(20))
                attributes.append(cdr3)
                attributes.append(str(quantity))
                attributes.append(str(len(cdr3)))
                prot = ProteinAnalysis(cdr3)
                attributes.append(f'{prot.molecular_weight():0.4f}')
                attributes.append(f'{prot.aromaticity():0.4f}')
                attributes.append(f'{prot.isoelectric_point():0.4f}')
                attributes.append(f'{aa_flex(cdr3)}')
                attributes.append(f'{prot.gravy():0.4f}')
                attributes.append(f'{prot.secondary_structure_fraction()[0]:0.4f}')
                attributes.append(f'{prot.secondary_structure_fraction()[1]:0.4f}')
                attributes.append(f'{prot.secondary_structure_fraction()[2]:0.4f}')

                # preferred to count the amino acids instead of
                # showing the percentage of them
                for num_of_fragment in prot.count_amino_acids().values():
                    attributes.append(str(num_of_fragment))

                # Here is the code to put percentage
                # for percent_of_fragment in prot.get_amino_acids_percent().values():
                    # attributes.append(f'{percent_of_fragment:0.4f}')

                groups = aa_groups(cdr3)
                for k, v in groups.items():
                    attributes.append(str(v))

                # in the final collumn we put the file
                # name where the data come from
                attributes.append(f'{output_file.name}')
                out.writerow(attributes)
                aa_total += 1

            except Exception:
                aa_error += 1
                # TODO: return ambiguous cdr3 sequences too?
                # print(f'\t{cdr3}\n')

    # TODO: return ambiguous cdr3 sequences too?
    # print(f"Total: {aa_total}\nError: {aa_error}\nPercentage: {aa_error/aa_total}")


def aa_flex(aa):
    total_flex = 0
    for residue in aa:
        for k, v in flex.items():
            if residue.upper() == k:
                total_flex += float(v[7])
    total_flex /= len(aa)
    return f'{total_flex:0.4f}'


flex = {
        'A': ['0.718', '0.717', '0.004', '0.556', '0.556', '0.009', '0.704', '0.704', '0.005', '0.847', '0.847', '0.009'],
        'C': ['0.668', '0.668', '0.010', '0.607', '0.607', '0.014', '0.671', '0.671', '0.010', '0.671', '0.670', '0.015'],
        'D': ['0.921', '0.921', '0.006', '0.726', '0.726', '0.011', '0.889', '0.889', '0.007', '1.055', '1.055', '0.014'],
        'E': ['0.963', '0.963', '0.005', '0.806', '0.805', '0.014', '0.912', '0.911', '0.005', '1.110', '1.110', '0.016'],
        'F': ['0.599', '0.599', '0.004', '0.465', '0.465', '0.010', '0.582', '0.582', '0.005', '0.653', '0.653', '0.011'],
        'G': ['0.843', '0.843', '0.005', '0.651', '0.651', '0.006', '0.811', '0.811', '0.009', '0.967', '0.967', '0.007'],
        'H': ['0.754', '0.754', '0.010', '0.598', '0.597', '0.009', '0.734', '0.734', '0.010', '0.894', '0.894', '0.014'],
        'I': ['0.632', '0.632', '0.004', '0.510', '0.510', '0.009', '0.617', '0.617', '0.004', '0.685', '0.686', '0.008'],
        'K': ['0.912', '0.912', '0.006', '0.863', '0.863', '0.007', '0.862', '0.862', '0.008', '1.016', '1.016', '0.009'],
        'L': ['0.681', '0.681', '0.003', '0.504', '0.504', '0.009', '0.650', '0.650', '0.007', '0.788', '0.788', '0.005'],
        'M': ['0.685', '0.685', '0.006', '0.575', '0.575', '0.014', '0.641', '0.641', '0.010', '0.740', '0.740', '0.013'],
        'N': ['0.851', '0.851', '0.008', '0.736', '0.735', '0.011', '0.848', '0.848', '0.009', '0.901', '0.901', '0.017'],
        'P': ['0.850', '0.850', '0.004', '0.753', '0.752', '0.006', '0.866', '0.866', '0.008', '0.857', '0.857', '0.009'],
        'Q': ['0.849', '0.849', '0.007', '0.730', '0.729', '0.007', '0.817', '0.817', '0.008', '1.007', '1.007', '0.015'],
        'R': ['0.814', '0.814', '0.006', '0.676', '0.676', '0.011', '0.807', '0.807', '0.006', '0.942', '0.942', '0.010'],
        'S': ['0.841', '0.840', '0.008', '0.698', '0.698', '0.014', '0.847', '0.846', '0.008', '0.915', '0.914', '0.010'],
        'T': ['0.758', '0.758', '0.004', '0.648', '0.648', '0.008', '0.742', '0.742', '0.006', '0.861', '0.862', '0.015'],
        'V': ['0.619', '0.619', '0.002', '0.503', '0.503', '0.006', '0.603', '0.603', '0.003', '0.707', '0.707', '0.009'],
        'W': ['0.627', '0.626', '0.011', '0.578', '0.577', '0.024', '0.609', '0.609', '0.013', '0.656', '0.656', '0.011'],
        'Y': ['0.615', '0.615', '0.004', '0.460', '0.461', '0.008', '0.567', '0.567', '0.005', '0.740', '0.741', '0.009']
        }


def handle_files(directory, pattern, type='VH'):
    for file in pathlib.Path(directory).rglob(pattern):
        pass


if __name__ == "__main__":
    import readline
    readline.set_completer_delims(' \t\n=')
    readline.parse_and_bind("tab: complete")
    try:

        start = timeit.default_timer()

        file = sys.argv[1]
        output = sys.argv[2]
        prefix = suffix = ''
        if len(sys.argv) > 3:
            change_name = sys.argv[3:]
            for modifier in change_name:
                if modifier.startswith('p='):
                    prefix = modifier[2:]
                elif modifier.startswith('s='):
                    suffix = modifier[2:]

        write_cdr3_attributes(file, output, prefix, suffix)

        print(f'Elapsed time: {timeit.default_timer() - start}\n\n')

    except KeyboardInterrupt:
        print('\n\nExiting ATTILA...\n')
