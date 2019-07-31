print('\n\n')
from Bio.SeqUtils.ProtParam import ProteinAnalysis

'''
# for fragment in y.count_amino_acids().items():
#   print(f'{fragment[1]},')

# for fragment in x.get_amino_acids_percent().items():
#   print(f'{fragment[1]:0.2f}')

# l = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

# cntX = []

# for f in l:
#     cntX.append('cnt'+f)
#     print(fr'nº{f}', end=',')

# print(y.molecular_weight())

# print(f'{y.aromaticity():0.4f},')
'''


x = ProteinAnalysis('FIVESK')
y = ProteinAnalysis('ATIDMGEYSSSSQWFDP')
prot = ProteinAnalysis('PKFPSNDAFDI')


sec_struc = prot.secondary_structure_fraction() # [helix, turn, sheet]

print(sec_struc)





print('\n\n')
