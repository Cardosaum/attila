import os
import time
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

start = time.time()

cdr3_VH_Regex = re.compile(r'((C)(\w+)(C)(..)(\w+)(WG.G))')

read = False

# regex for lines that starts with '#'
hashRegex = re.compile(r'(([\d\w]+)(:)(\d+)(:)(\d+)(-)([\d\w]+)(:)(\d)(:)(\d+)(:)(\d+)(:)(\d+)([|])(FRAME:(\d))([|])(\[(\d+)-(\d+)\])([|])(\d+)([|])(\d+)(\.)?(\d+))')

# regex for lines that starts with '*'
asterixRegex = re.compile(r'((C)(\w+)(C)(.{2})(\w+)(WG.G))')

# regex for lines that starts with '>'
greaterRegex = re.compile(r'((\w+)(:)(\w+)(:)(\d+)(-)(\w+)(:)(\d)(:)(\d+)(:)(\d+)(:)(\d+)([|])(FRAME:(\d+))(\s+)?(\w+)?)')

# regex for lines that does not start starts with an special character
seqRegex = re.compile(r'((.+)(C)(.+)(C)(.{2})(.+)(WG.G)(.+)?)')

seqPatterns = {}

# lisFiles = ['/home/matheus/mcs/wo/R0/Renato_zika_acido_R0_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/rafaCD20_Vh_R0_R1aafreq.txt', '/home/matheus/mcs/wo/R0/Renato__zika_R0_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR2_S2_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR4_S8_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR2_S6_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR0_S1_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR4_S4_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR3_S7_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL29VHR3_S3_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/Thais_29_66/VCL66VHR0_S5_L001_R1_001aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/gal66R0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/gal29R0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/cd20rafaR0aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/R4pep_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/R4ac_VH_R1aafreq.txt', '/home/matheus/mcs/wo/R0/analiseR0/zika/zikaR0aafreq.txt']

lisFiles = ['/home/matheus/mcs/wo/R0/Renato_zika_acido_R0_VH_R1aafreq.txt']

for ffile in lisFiles:

    outputFile = ffile.split('/')[-1].replace('aafreq.txt', '.csv')

    with open(ffile, encoding='ISO-8859-1') as file, open(f'/home/matheus/Documentos/test/{outputFile}', 'w') as out:

        # creat the header for outputFile
        out.write('id,cdr3,length,cntA,nºA,cntC,nºC,cntD,nºD,cntE,nºE,cntF,nºF,cntG,nºG,cntH,nºH,cntI,nºI,cntK,nºK,cntL,nºL,cntM,nºM,cntN,nºN,cntP,nºP,cntQ,nºQ,cntR,nºR,cntS,nºS,cntT,nºT,cntV,nºV,cntW,nºW,cntY,nºY,pctA,N_pctA,pctC,N_pctC,pctD,N_pctD,pctE,N_pctE,pctF,N_pctF,pctG,N_pctG,pctH,N_pctH,pctI,N_pctI,pctK,N_pctK,pctL,N_pctL,pctM,N_pctM,pctN,N_pctN,pctP,N_pctP,pctQ,N_pctQ,pctR,N_pctR,pctS,N_pctS,pctT,N_pctT,pctV,N_pctV,pctW,N_pctW,pctY,N_pctY,\n')

        for line in file:
          if line[0] == '#':
            try:
                parse = hashRegex.search(line)
                sequenceID = parse.groups()[0]
                # print(sequenceID)
            except:
              print(f'ERROR: {ffile}, line: {line}')
          elif line[0] == '*':
              parse = asterixRegex.search(line)
          elif line[0] == '>':
              try:
                  parse = greaterRegex.search(line)
                  parse = parse[0].replace('\n', '')
                  # print(parse[0])
                  out.write(f'{parse},')
              except:
                  print(f'ERROR: {ffile}, line: {line}')
              read = True
          elif read:
              parse = seqRegex.search(line)
              parseLen = str(len(parse[7]))
              out.write(f'{parse[7]},')
              out.write(f'{parseLen},')
              prot = ProteinAnalysis(parse[7])
              for fragment in prot.count_amino_acids().items():
                out.write(f'{fragment[0]},{fragment[1]},')
              for fragment in prot.get_amino_acids_percent().items():
                out.write(f'{fragment[0]},{fragment[1]:0.2f},')
              out.write('\n')
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

end = time.time()

print('Elapsed time: ' + str(end - start))
