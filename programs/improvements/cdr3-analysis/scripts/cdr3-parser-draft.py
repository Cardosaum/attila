# print('\n\n')
# from Bio.SeqUtils.ProtParam import ProteinAnalysis

# '''
# # for fragment in y.count_amino_acids().items():
# #   print(f'{fragment[1]},')

# # for fragment in x.get_amino_acids_percent().items():
# #   print(f'{fragment[1]:0.2f}')

# # l = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

# # cntX = []

# # for f in l:
# #     cntX.append('cnt'+f)
# #     print(fr'nº{f}', end=',')

# # print(y.molecular_weight())

# # print(f'{y.aromaticity():0.4f},')

# sec_struc = prot.secondary_structure_fraction() # [helix, turn, sheet]

# print(sec_struc)

# epsilon_prot = prot.molar_extinction_coefficient()  # [Reduced, Oxidized]

# print(epsilon_prot)
# '''


# x = ProteinAnalysis('FIVESK')
# w = ProteinAnalysis('MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV')
# y = ProteinAnalysis('ATIDMGEYSSSSQWFDP')
# prot = ProteinAnalysis('PKFPSNDAFDI')

# d = ProteinAnalysis('GPPHPNYYYHMDV')

# print(f'{d.molecular_weight():0.4f}')

# print('\n\n')

lis = ['cdr3', 'length', 'MW', 'AV', 'II', 'IP', 'SSF_Helix', 'SSF_Turn', 'SSF_Sheet', 'MEC_Reduced', 'MEC_Oxidized', 'nºA', 'nºC', 'nºD', 'nºE', 'nºF', 'nºG', 'nºH', 'nºI', 'nºK', 'nºL', 'nºM', 'nºN', 'nºP', 'nºQ', 'nºR', 'nºS', 'nºT', 'nºV', 'nºW', 'nºY', '%A', '%C', '%D', '%E', '%F', '%G', '%H', '%I', '%K', '%L', '%M', '%N', '%P', '%Q', '%R', '%S', '%T', '%V', '%W', '%Y']

out = ''
for item in lis:
    out += f'{item},'

print(out)
