import os
import time
import re


cdr3_VH_Regex = re.compile(r'((C)(\w+)(C)(..)(\w+)(WG.G))')

read = False

hashRegex = re.compile(r'(([\d\w]+)(:)(\d+)(:)(\d+)(-)([\d\w]+)(:)(\d)(:)(\d+)(:)(\d+)(:)(\d+)([|])(FRAME:(\d))([|])(\[(\d+)-(\d+)\])([|])(\d+)([|])(\d+)(\.)(\d+))')

asterixRegex = re.compile(r'((C)(\w+)(C)(.{2})(\w+)(WG.G))')

greaterRegex = re.compile(r'((\w+)(:)(\w+)(:)(\d+)(-)(\w+)(:)(\d)(:)(\d+)(:)(\d+)(:)(\d+)([|])(FRAME:(\d+))(\s+)(\w+))')

seqRegex = re.compile(r'((.+)(C)(.+)(C)(.{2})(.+)(WG.G)(.+))')

with open('/home/matheus/mcs/wo/R0/Renato_zika_acido_R0_VH_R1aafreq.txt', encoding='ISO-8859-1') as file:
    for line in file:
      if line[0] == '#':
          parse = hashRegex.search(line)
          # TODO: find what is the usefull data here
      elif line[0] == '*':
          parse = asterixRegex.search(line)
      elif line[0] == '>':
          parse = greaterRegex.search(line)
          read = True
      elif read:
          parse = seqRegex.search(line)
          print(parse[7])
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
