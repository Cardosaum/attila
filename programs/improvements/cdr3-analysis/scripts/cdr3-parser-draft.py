import os
import re

dirToWatch = '/home/matheus/mcs/wo/R0/'

fileRegex = re.compile(r'((.+)aafreq.txt)')

filesInDir = []

for dirPath, dirName, filenames in os.walk(dirToWatch):
    for file in filenames:
      # print(os.path.join(dirPath, file))
      match = re.search(fileRegex, file)
      if match and file[0] != '.':
        print(os.path.join(dirPath, file))
        filesInDir.append(os.path.join(dirPath, file))
print(filesInDir)
