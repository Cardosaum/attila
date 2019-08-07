import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# regex for lines that does not start starts with an special character
seqRegex = re.compile(r'((.+)(C)(.+)(C)(.{2})(.+)(WG.G)(.+)?)')

def cdr3Analysis(cdr3Sequence):
    parse = seqRegex.search(cdr3Sequence)
    print(parse)
