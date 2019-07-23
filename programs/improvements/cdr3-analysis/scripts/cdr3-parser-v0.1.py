import os
import time
import re

with open('/home/matheus/mcs/study/code/bioinfo/attila/programs/improvements/cdr3-analysis/files/R0/test01.txt') as file:
    for line in file:
        if line[0] == '*':
            print(line)
            for char in range(len(line)):
                if line[char] == 'C':
                    print(line[char], 'at index: ', char)



                    print('So, CDR3 is: {}'.format(line[char+2]))

