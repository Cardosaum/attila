#!/usr/bin/env python3
'''

cdr3Parser.py - Analyse an "aafreq.txt" input file and return all features of CDR3 sequences

Usage:
       python3 cdr3Parser.py <INPUT_FILE>

Copyright:  (c)  2019  Matheus Cardoso   <https://github.com/cardosaum>

License:  Apache 2.0  <https://www.apache.org/licenses/LICENSE-2.0>

'''


import os
import time
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# To show elapsed time
start = time.time()

cdr3_VH_Regex = re.compile(r'((C)(\w+)(C)(..)(\w+)(WG.G))')

# Needed in order to verify if the line is or isn't the actual full sequence,
# like: "YAAQPAMAEVQLLESGGGVVQPGRSLRLSCVASGFNFSPYWMHWVRQAPGKGLVWVSHIHSDGTSTSYADSVKGRFTISRDNAKNTLYLEMNSLRPEDTAVYYCVTFIVESKWGQGTLITVSSASTKGPS", for example
read = False

# regex for lines that starts with '#'
# this regex is for analysing sequences' IDs
# The expected sequences' IDs are, for example:
# "#M04816:19:000000000-D4997:1:1102:14947:28739|FRAME:3|[36-122]|84717|1227616653092.06"
# You can use Regex editors such as <https://regex101.com/> for a better understanding
# Just copy and paste this regex and sequence ID to see how it's working
hashRegex = re.compile(r'(([\d\w]+)(:)(\d+)(:)(\d+)(-)([\d\w]+)(:)(\d)(:)(\d+)(:)(\d+)(:)(\d+)([|])(FRAME:(\d))([|])(\[(\d+)-(\d+)\])([|])(\d+)([|])(\d+)(\.)?(\d+))')

# regex for lines that starts with '*'
asterixRegex = re.compile(r'((C)(\w+)(C)(.{2})(\w+)(WG.G))')

# regex for lines that starts with '>'
greaterRegex = re.compile(r'((\w+)(:)(\w+)(:)(\d+)(-)(\w+)(:)(\d)(:)(\d+)(:)(\d+)(:)(\d+)([|])(FRAME:(\d+))(\s+)?(\w+)?)')

# regex for lines that does not start starts with an special character
seqRegex = re.compile(r'((.+)(C)(.+)(C)(.{2})(.+)(WG.G)(.+)?)')

# TODO: Is this really necessary?
seqPatterns = {}

# TODO: Remove this list when finish to write the script
# lisFiles = ['/home/matheus/mcs/wo/R0/Renato_zika_acido_R0_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/rafaCD20_Vh_R0_R1aafreq.txt', '/home/matheus/mcs/wo/R0/Renato__zika_R0_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR2_S2_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR4_S8_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR2_S6_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR0_S1_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR4_S4_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR3_S7_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR3_S3_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR0_S5_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/gal66R0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/gal29R0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/cd20rafaR0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/R4pep_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/R4ac_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/zikaR0aafreq.txt']

# TODO: Remove this list when finish to write the script
lisFiles = ['/home/matheus/mcs/wo/R0/Renato_zika_acido_R0_VH_R1aafreq.txt']

# TODO: Loop for all items in "lisFiles" and analyze each one separately
for ffile in lisFiles:

    # Get the "Base Name" of the input file, replacing "aafreq.txt" for ".csv"
    outputFile = ffile.split('/')[-1].replace('aafreq.txt', '.csv')

    # This encoding is needed in order to be compatible with different files
    # In previous analyses was found that some "aafreq.txt" files had
    # strange characters (Error in sequencing, maybe?)
    with open(ffile, encoding='ISO-8859-1') as file, open(f'/home/matheus/Documentos/test/{outputFile}', 'w') as out:

        # creat the header for outputFile
        out.write(r'id,cdr3,length,MW,AV,II,IP,SSF_Helix,SSF_Turn,SSF_Sheet,MEC_Reduced,MEC_Oxidized,nºA,nºC,nºD,nºE,nºF,nºG,nºH,nºI,nºK,nºL,nºM,nºN,nºP,nºQ,nºR,nºS,nºT,nºV,nºW,nºY,%A,%C,%D,%E,%F,%G,%H,%I,%K,%L,%M,%N,%P,%Q,%R,%S,%T,%V,%W,%Y' + '\n')

        # Analyze each line in order to know what is the content
        for line in file:
          # Expected pattern is something similar to this:
          # "#M04816:19:000000000-D4997:1:1102:14947:28739|FRAME:3|[36-122]|84717|1227616653092.06"
          if line[0] == '#':

            # Try to get all sequenceID
            # TODO: Get something more, or one variable for each sub-pattern matched?
            try:
                parse = hashRegex.search(line)
                sequenceID = parse.groups()[0]
                # TODO: Remove when script is finished
                # print(sequenceID)
            except:
              print(f'ERROR: {ffile}, line: {line}')

          # Expected pattern is something similar to this:
          # "*CVASGFNFSPYWMHWVRQAPGKGLVWVSHIHSDGTSTSYADSVKGRFTISRDNAKNTLYLEMNSLRPEDTAVYYCVTFIVESKWGQG"
          elif line[0] == '*':
              parse = asterixRegex.search(line)

          # Expected pattern is something similar to this:
          # ">M04816:19:000000000-D4997:1:1101:19065:1972|FRAME:3"
          elif line[0] == '>':

              # Try to get informations about sequenceID
              # and write this informations to the output file
              # TODO: Get something more, or one variable for each sub-pattern matched?
              try:
                  parse = greaterRegex.search(line)
                  parse = parse[0].replace('\n', '')
                  # print(parse)
                  out.write(f'{parse},')
              except:
                  print(f'ERROR: {ffile}, line: {line}')
              # This line is needed to inform that the next line will be the protein sequence of the previous analysed sequence ID
              read = True

          # If the previous line in <INPUT_FILE> was a line starting with ">", then, the next line will be the respective protein sequence
          elif read:
              parse = seqRegex.search(line)

              # Here, parse[7] correspond to the CDR3 sequence
              out.write(f'{parse[7]},')

              # We're getting the length of CDR3 sequence
              parseLen = str(len(parse[7]))
              out.write(f'{parseLen},')

              # We'll start to analyze the CDR3 sequence
              # For more information about each function,
              # Check out <http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam-module.html>
              prot = ProteinAnalysis(parse[7])

              # Get Molecular weight ot CDR3 sequence
              out.write(f'{prot.molecular_weight()},')

              # Get aromaticity value
              out.write(f'{prot.aromaticity():0.4f},')

              # Get instability index
              out.write(f'{prot.instability_index():0.4f},')

              # Get isoelectirc point
              out.write(f'{prot.isoelectric_point():0.4f},')

              # Get Secondary Structure Information
              sec_struc = prot.secondary_structure_fraction() # [helix, turn, sheet]

              sec_struc_helix = sec_struc[0]
              sec_struc_turn = sec_struc[1]
              sec_struc_sheet = sec_struc[2]

              out.write(f'{sec_struc_helix:0.4f},')
              out.write(f'{sec_struc_turn:0.4f},')
              out.write(f'{sec_struc_sheet:0.4f},')

              # Get Molar Extinction Coefficient
              # TODO: Understand what this classification means
              epsilon_prot = prot.molar_extinction_coefficient()  # [Reduced, Oxidized]
              epsilon_prot_reduced = epsilon_prot[0]  # with reduced cysteines
              epsilon_prot_oxidized = epsilon_prot[1]  # with disulfid bridges

              out.write(f'{epsilon_prot_reduced},')
              out.write(f'{epsilon_prot_oxidized},')

              # Get the aminoacids count
              for fragment in prot.count_amino_acids().items():
                out.write(f'{fragment[1]},')

              # Get their percentage
              for fragment in prot.get_amino_acids_percent().items():
                out.write(f'{fragment[1]:0.4f},')

              out.write('\n')
              # Here, we are ensuring that the loop will only execute if reach another line beginning with ">"
              read = False


# TODO:
''' Find which properties should be analised
-Modules from EMBOSS

- inforesidue
- pepstats
- pepinfo
- charge
- hmoment
- iep
- octanol
- pepwindow
- pepwindowall
- garnier?
'''

# Here we get when the program finished
end = time.time()

# Print the elapsed time, in seconds
print('Elapsed time: ' + str(end - start))
