from Bio.SeqUtils.ProtParam import ProteinAnalysis

x = ProteinAnalysis('FIVESK')
y = ProteinAnalysis('XAQPAXAQVQLXESGXGLVQPGGSLRLSCVASGFDFSRYWMHWVRQAPGKGLEWVSHIHSDGIPTAYADSVRGRFTISRDISKNTLYLQMNNLRPEDTAVYYCVTFIVESKWGQGTLATVSSASTXGPSYSXHAXX')

# for fragment in y.count_amino_acids().items():
#   print(f'{fragment[0]},{fragment[1]}', end=',')

for fragment in x.get_amino_acids_percent().items():
  print(f'{fragment[0]},{fragment[1]:0.2f}', end=',')

l = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

# cntX = []

# for f in l:
#   cntX.append('cnt'+f)
#   print(f'pct{f},N_pct{f}', end=',')
