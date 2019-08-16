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



all_cdr3 = extract_cdr3('/home/mcsouza/mcs-uracila/files-to-work-with/R0/analiseR0/zika/R4pep_VH_R1aafreq.txt')

# print(all_cdr3)

# Write to a file the attributes of CDR3 sequence
# In order, the attributes written are:
# TODO: Place HEADER here
def get_cdr3_attributes(file):
  # TODO: Write HEADER to file
  output_file = f'{os.path.splitext(file)[0]}.csv'
  with open(output_file, 'w') as out:
    out.write(r'cdr3;quantity;length;MW;AV;IP;flex;gravy;SSF_Helix;SSF_Turn;SSF_Sheet;nºA;nºC;nºD;nºE;nºF;nºG;nºH;nºI;nºK;nºL;nºM;nºN;nºP;nºQ;nºR;nºS;nºT;nºV;nºW;nºY;%A;%C;%D;%E;%F;%G;%H;%I;%K;%L;%M;%N;%P;%Q;%R;%S;%T;%V;%W;%Y' + '\n')
    for cdr3, quantity in all_cdr3.items():
      attributes = []
      # print(f'CDR3:\t{cdr3}'.ljust(60)+f'Quantity: {quantity}'.rjust(20))
      attributes.append(cdr3)
      attributes.append(str(quantity))
      attributes.append(str(len(cdr3)))
      prot = ProteinAnalysis(cdr3)
      attributes.append(f'{prot.molecular_weight():0.4f}')
      attributes.append(f'{prot.aromaticity():0.4f}')
      attributes.append(f'{prot.isoelectric_point():0.4f}')
      attributes.append(f'{prot.flexibility()}')
      attributes.append(f'{prot.gravy():0.4f}')
      attributes.append(f'{prot.secondary_structure_fraction()[0]:0.4f}')
      attributes.append(f'{prot.secondary_structure_fraction()[1]:0.4f}')
      attributes.append(f'{prot.secondary_structure_fraction()[2]:0.4f}')
      for num_of_fragment in prot.count_amino_acids().values():
        attributes.append(str(num_of_fragment))
      for percent_of_fragment in prot.get_amino_acids_percent().values():
        attributes.append(f'{percent_of_fragment:0.4f}')
      # TODO: ASA (accessibility) - (?)
      # TODO: Write all attributes to a file
      out.write(f'{";".join(attributes)}\n')

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
	print(len(sequence))
	r = 0
	for k, v in result.items():
		if not k == 'invalid':
				r += v
	print(r)
	print(not_listed)
	return result





get_cdr3_attributes('/home/mcsouza/mcs-uracila/files-to-work-with/R0/analiseR0/zika/R4pep_VH_R1aafreq.txt')

print(f'\n\nElapsed time: {timer() - start}')

print(aa_groups('XÿGMÿGVAAQPAMAQVQLQESGGGLVQPGGSLRLSCVASGFDFSRYWMHWVRQAPGKGLEWVSHIHSDGIPTAYADSVRGRFTISRDISKNTLYLQMNNLRPEDTAVYYCVTFIVESKWGQGTLVTVSSAXTKGPS'))
