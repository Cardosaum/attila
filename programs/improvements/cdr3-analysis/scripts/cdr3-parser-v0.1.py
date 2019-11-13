#!/usr/bin/env python3
'''

cdr3Parser.py - Analyse an "aafreq.txt" input file and return all features of CDR3 sequences

Usage:
       python3 cdr3Parser.py <INPUT_FILE>

Copyright:  (c)  2019  Matheus Cardoso   <https://github.com/cardosaum>

License:  Apache 2.0  <https://www.apache.org/licenses/LICENSE-2.0>

'''


import os
import sys
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import defaultdict
from timeit import default_timer as timer

# To show elapsed time
start = timer()

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

# Creat a dictionary for count CDR3 sequences that are repeated
countCDR3 = defaultdict(int)

# TODO: Remove this list when finish to write the script

dirFiles = sys.argv[1:][0].strip()
tmpFiles = os.listdir(dirFiles)
lisFiles = []
for file in tmpFiles:
  if file.startswith('R'):
    lisFiles.append(file)
lisFiles.sort()



# TODO: Loop for all items in "lisFiles" and analyze each one separately
for ffile in lisFiles:
    # Get the "Base Name" of the input file, replacing "aafreq.txt" for ".csv"
    # outputFile = ffile.split('/')[-1].replace('aafreq.txt', '.csv')
    ffile = os.path.join(dirFiles, ffile)
    print(ffile)
    outputFile = os.path.basename(ffile).replace('aafreq.txt', 'aafreq.csv')
    outputPath = os.path.dirname(ffile)
    outputAbsPath = os.path.join(outputPath, outputFile)

    # This encoding is needed in order to be compatible with different files
    # In previous analyses was found that some "aafreq.txt" files had
    # strange characters (Error in sequencing, maybe?)
    with open(ffile) as file, open(outputAbsPath, 'w') as out:

        # Legend for the header in output file:
        #
        # NOTE: All features above were evaluated using Biopython mudule,
        #       For more information about each feature, read the docs.
        #
        #       ProteinAnalysis - Biopython:
        #       <http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam-module.html>
        #
        #         (all subsequent characteristics correspont to cdr3 sequence)
        #
        # cdr3: Is the CDR3 sequence of the protein fragment.
        #       We considered CDR3 the sequence with the folowing pattern:
        #
        #       <random_sequence>C<random_sequence>Cxx<cdr3_sequence>WGxG
        #
        #       Where "x" stands for any random aminoacid
        #
        # nºAPP: Number of appearances - Number of different protein sequences
        #        That have same CDR3 sequence;
        #
        #        Add "1" for each unique sequence found in analysis
        #
        # length: Length of the given CDR3 sequence
        #
        # MW: Molecular Weight
        #
        # AV: Aromaticity Value
        #
        # II: Instability Index
        #
        # IP: Isoelectric Point
        #
        # SSF_Helix: Secondary Structure Fragment - Helix
        #
        # SSF_Turn: Secondary Structure Fragment - Turn
        #
        # SSF_Sheet: Secondary Structure Fragment - Sheet
        #
        # MEC_Reduced: Molar Extinction Coefficient - Reduced, (with reduced cysteines)
        #
        # MEC_Oxidized: Molar Extinction Coefficient - Oxidized, (with dissulfid bridges)
        #
        # nºX: Number of aminoacids X in CDR3 sequence
        #      So, if the sequence are "AAACCD"
        #      nºA = 3, nºC = 2, nºD = 1, nºE = 0, nºF = 0, and so on
        #
        # %X: Percentage of aminoacid X in CDR3 sequence
        #     The logic is similar to "nºX"

        # creat the header for outputFile
        out.write(r'cdr3,length,MW,AV,II,IP,SSF_Helix,SSF_Turn,SSF_Sheet,MEC_Reduced,MEC_Oxidized,nºA,nºC,nºD,nºE,nºF,nºG,nºH,nºI,nºK,nºL,nºM,nºN,nºP,nºQ,nºR,nºS,nºT,nºV,nºW,nºY,%A,%C,%D,%E,%F,%G,%H,%I,%K,%L,%M,%N,%P,%Q,%R,%S,%T,%V,%W,%Y'
+ '\n')

        # Analyze each line in order to know what is the content
        for line in file:

          # List of all attributes that will be written on outputfile
          listOfAttributes = []

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
              # try:
              #     parse = greaterRegex.search(line)
              #     parse = parse[0].replace('\n', '')
              #     # print(parse)
              #     out.write(f'{parse},')
              # except:
              #     print(f'ERROR: {ffile}, line: {line}')
              # This line is needed to inform that the next line will be the protein sequence of the previous analysed sequence ID
              read = True

          # If the previous line in <INPUT_FILE> was a line starting with ">", then, the next line will be the respective protein sequence
          elif read:

              parse = seqRegex.search(line)

              cdr3_sequence = parse[7]
              # Analyzes the sequence and store the Regex Object
              try:
                  prot = ProteinAnalysis(parse[7])
              except:
                print(f'ERRO: {ffile}, line: {line}')

              if cdr3_sequence in countCDR3.keys():
                  # If CDR3 sequence already have been analyzed, only add "1"
                  # in "countCDR3"

                  # Add "1" for each unique CDR3 sequence
                  countCDR3[cdr3_sequence] += 1
                  read = False

              else:
                  # Add "1" for each unique CDR3 sequence
                  countCDR3[cdr3_sequence] += 1

                  # # Here, parse[7] correspond to the CDR3 sequence
                  # try:
                  #     listOfAttributes.append(f'{parse[7]}')
                  # except:
                  #   print(f'ERRO: {ffile}, line: {line}')

                  # Add unique CDR3 sequence to outputfile
                  listOfAttributes.append(cdr3_sequence)

                  # We're getting the length of CDR3 sequence
                  parseLen = str(len(parse[7]))
                  try:
                      listOfAttributes.append(f'{parseLen}')
                  except:
                    print(f'ERRO: {ffile}, line: {line}')

                  # We'll start to analyze the CDR3 sequence
                  # For more information about each function,
                  # Check out <http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam-module.html>


                  # Get Molecular weight ot CDR3 sequence
                  try:
                      listOfAttributes.append(f'{prot.molecular_weight():0.4f}')
                  except:
                    print(f'ERRO: {ffile}, line: {line}')

                  # Get aromaticity value
                  try:
                      listOfAttributes.append(f'{prot.aromaticity():0.4f}')
                  except:
                    print(f'ERRO: {ffile}, line: {line}')

                  # Get instability index
                  try:
                      listOfAttributes.append(f'{prot.instability_index():0.4f}')
                  except:
                    print(f'ERRO: {ffile}, line: {line}')

                  # Get isoelectirc point
                  try:
                      listOfAttributes.append(f'{prot.isoelectric_point():0.4f}')
                  except:
                    print(f'ERRO: {ffile}, line: {line}')

                  # Get Secondary Structure Information
                  try:
                      sec_struc = prot.secondary_structure_fraction() # [helix, turn, sheet]

                      sec_struc_helix = sec_struc[0]
                      sec_struc_turn = sec_struc[1]
                      sec_struc_sheet = sec_struc[2]
                  except:
                    print(f'ERRO: {ffile}, line: {line}')

                  try:
                      listOfAttributes.append(f'{sec_struc_helix:0.4f}')
                      listOfAttributes.append(f'{sec_struc_turn:0.4f}')
                      listOfAttributes.append(f'{sec_struc_sheet:0.4f}')
                  except:
                    print(f'ERRO: {ffile}, line: {line}')

                  # Get Molar Extinction Coefficient
                  # TODO: Understand what this classification means
                  try:
                      epsilon_prot = prot.molar_extinction_coefficient()  # [Reduced, Oxidized]
                      epsilon_prot_reduced = epsilon_prot[0]  # with reduced cysteines
                      epsilon_prot_oxidized = epsilon_prot[1]  # with disulfid bridges

                      listOfAttributes.append(f'{epsilon_prot_reduced}')
                      listOfAttributes.append(f'{epsilon_prot_oxidized}')
                  except:
                      print(f'ERRO: {ffile}, line: {line}')


                  # Get the aminoacids count
                  try:
                      for fragment in prot.count_amino_acids().items():
                          listOfAttributes.append(f'{fragment[1]}')
                  except:
                    print(f'ERRO: {ffile}, line: {line}')

                  # Get their percentage
                  try:
                      for fragment in prot.get_amino_acids_percent().items():
                          listOfAttributes.append(f'{fragment[1]:0.4f}')
                  except:
                    print(f'ERRO: {ffile}, line: {line}')

                  # Write all attributes of the current CDR3 sequence to
                  # output file
                  writeAttributes = ','.join(listOfAttributes)
                  out.write(f'{writeAttributes}\n')
                  # print(f'{writeAttributes}\n')

                  # Here, we are ensuring that the loop will only execute if
                  # reach another line beginning with ">"
                  read = False

# print(countCDR3)

# Here we get when the program finished
end = timer()

# Print the elapsed time, in seconds
print(f'Elapsed time: {end - start}')
