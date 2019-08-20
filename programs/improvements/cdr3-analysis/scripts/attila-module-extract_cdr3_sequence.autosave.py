#!/usr/bin/env python3
'''

attila-module-extract_cdr3_sequence.py - Analyze an "aafreq.txt" input file and return all features of CDR3 sequences

Usage:
       python3 cdr3Parser.py <INPUT_FILE>

Copyright:  (c)  2019  Matheus Cardoso   <https://github.com/cardosaum>

License:  Apache 2.0  <https://www.apache.org/licenses/LICENSE-2.0>

'''

import re
import os
from timeit import default_timer as timer
from Bio.SeqUtils.ProtParam import ProteinAnalysis

start = timer()

cdr3_regex = re.compile(r'((.+)(C)(.+)(C)(.{2})(.+)(WG.G)(.+)?)')

def extract_cdr3(file):
  dic_cdr3 = {}
  with open(file, encoding="ISO-8859-1") as f:
    read = False
    for line in f:
      if line.startswith('>'):
        read = True
      elif read:
        parse = cdr3_regex.search(line)
        cdr3 = parse[7]
        if not cdr3 in dic_cdr3:
          dic_cdr3.setdefault(cdr3, 1)
        else:
          dic_cdr3[cdr3] += 1
        read = False
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
	# print(len(sequence))
	r = 0
	for k, v in result.items():
		if not k == 'invalid':
				r += v
	# print(r)
	# print(not_listed)
	return result

# Write to a file the attributes of CDR3 sequence
# In order, the attributes written are:
# TODO: Place HEADER here
def write_cdr3_attributes(file):
  # TODO: Write HEADER to file
  output_file = f'{os.path.splitext(file)[0]}.csv'
  with open(output_file, 'w') as out:
    out.write(r'cdr3;quantity;length;MW;AV;IP;flex;gravy;SSF_Helix;SSF_Turn;SSF_Sheet;%A;%C;%D;%E;%F;%G;%H;%I;%K;%L;%M;%N;%P;%Q;%R;%S;%T;%V;%W;%Y;aliphatic;aromatic;neutral;positive;negative;invalid' + '\n')
    for cdr3, quantity in extract_cdr3(file).items():
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
      # for num_of_fragment in prot.count_amino_acids().values():
      #   attributes.append(str(num_of_fragment))
      for percent_of_fragment in prot.get_amino_acids_percent().values():
        attributes.append(f'{percent_of_fragment:0.4f}')
      groups = aa_groups(cdr3)
      for k, v in groups.items():
        attributes.append(str(v))
      # TODO: ASA (accessibility) - (?)
      # TODO: Write all attributes to a file
      out.write(f'{";".join(attributes)}\n')


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


lis_files = ['/home/matheus/mcs/wo/R0/._Renato__zika_R0_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/._Renato_zika_acido_R0_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/._rafaCD20_Vh_R0_R1aafreq.txt', '/home/matheus/mcs/wo/R0/Renato__zika_R0_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/Renato_zika_acido_R0_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/rafaCD20_Vh_R0_R1aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/._VCL29VHR0_S1_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/._VCL29VHR3_S3_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/._VCL66VHR0_S5_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/._VCL66VHR2_S6_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/._VCL66VHR3_S7_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/._VCL66VHR4_S8_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR0_S1_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR2_S2_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR3_S3_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR4_S4_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR0_S5_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR2_S6_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR3_S7_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR4_S8_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/._cd20rafaR0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/._gal29R0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/._gal66R0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/cd20rafaR0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/gal29R0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/gal66R0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/._R4ac_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/._R4pep_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/._zikaR0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/R4ac_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/R4pep_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/zikaR0aafreq.txt']


# write_cdr3_attributes()
for file in lis_files:
    print(file.replace('/home/matheus', os.path.expanduser('~')))

print(f'\n\nElapsed time: {timer() - start}')
