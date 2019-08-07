#!/usr/bin/env python3
'''

attila-module-extract_cdr3_sequence.py - Analyze an "aafreq.txt" input file and return all features of CDR3 sequences

Usage:
       python3 cdr3Parser.py <INPUT_FILE>

Copyright:  (c)  2019  Matheus Cardoso   <https://github.com/cardosaum>

License:  Apache 2.0  <https://www.apache.org/licenses/LICENSE-2.0>

'''

import re
import collections
import os
from timeit import default_timer as timer
from Bio.SeqUtils.ProtParam import ProteinAnalysis

start = timer()

cdr3_regex = re.compile(r'((.+)(C)(.+)(C)(.{2})(.+)(WG.G)(.+)?)')


def extract_cdr3(file):
  dic_cdr3 = {}
  with open(file) as f:
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



all_cdr3 = extract_cdr3('/home/matheus/mcs/wo/R0/analiseR0/zika/R4ac_VH_R1aafreq.txt')

# print(all_cdr3)

# Write to a file the attributes of CDR3 sequence
# In order, the attributes written are:
# TODO: Place HEADER here
def get_cdr3_attributes(file):
  # TODO: Write HEADER to file
  for cdr3, quantity in all_cdr3.items():
    attributes = []
    print(f'CDR3:\t{cdr3}'.ljust(60)+f'Quantity: {quantity}'.rjust(20))
    attributes.append(cdr3)
    attributes.append(quantity)
    attributes.append(len(cdr3))
    prot = ProteinAnalysis(cdr3)
    attributes.append(f'{prot.molecular_weight():0.4f}')
    attributes.append(f'{prot.aromaticity():0.4f}')
    attributes.append(f'{prot.isoelectric_point():0.4f}')
    attributes.append(f'{prot.secondary_structure_fraction()[0]:0.4f}')
    attributes.append(f'{prot.secondary_structure_fraction()[1]:0.4f}')
    attributes.append(f'{prot.secondary_structure_fraction()[2]:0.4f}')
    for num_of_fragment in prot.count_amino_acids().values():
      attributes.append(num_of_fragment)
    for percent_of_fragment in prot.get_amino_acids_percent().values():
      attributes.append(f'{percent_of_fragment:0.4f}')
    print(attributes)


print(get_cdr3_attributes('/home/matheus/mcs/wo/R0/analiseR0/zika/R4ac_VH_R1aafreq.txt'))

print(f'\n\nElapsed time: {timer() - start}')
