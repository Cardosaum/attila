import pprint

attila_commands = {
					'commands:':
						{
							'CTRL-C':'quit ATTILA; abort analysis',
							'TAB':'autocomplete a path'
						},
					'configuration parameters:':
						{
							'Configuration files exist (y or n)': 'type "y" if you already have configuration files type "n" or press ENTER key if you prefer to let ATTILA create the configuration files',
							'Path of the configuration file of VH libraries': 'location of the configuration file of VH libraries',
							'Path of the configuration file of VL libraries': 'location of the configuration file of VL libraries',
							'Project Name': 'name of the directory that will be created by ATTILA to save output files',
							'Directory to save project': 'the directory where the project will be saved',
							'Reads are paired-end (y or n)': 'type y or press ENTER key for yes; type n for no',
							'Minimum read length': 'default value is 300 pb; type y to change default; type n or press ENTER key to use default value if you choose to change default value, the new read length must be an integer number',
							'Minimum base quality': 'default value is 20; type y to change default; type n or press ENTER key to use default value if you choose to change default value, the new base quality must be an integer number',
							'Number of candidates to rank': 'number of candidate clones that ATTILA will try to find in VH and VL libraries the number must be an integer',
						},
					'Parameters for paired-end reads:':
						{
							'Path of fastq file of VH R0 reads r1': 'location of the fastq file containing reads r1 from initial VH library',
							'Path of fastq file of VH R0 reads r2': 'location of the fastq file containing reads r2 from initial VH library',
							'Path of fastq file of VH RN reads r1': 'location of the fastq file containing reads r1 from final VH library',
							'Path of fastq file of VH RN reads r2': 'location of the fastq file containing reads r2 from final VH library',
							'Path of fastq file of VL R0 reads r1': 'location of the fastq file containing reads r1 from initial VL library',
							'Path of fastq file of VL R0 reads r2': 'location of the fastq file containing reads r2 from initial VL library',
							'Path of fastq file of VL RN reads r1': 'location of the fastq file containing reads r1 from final VL library',
							'Path of fastq file of VL RN reads r2': 'location of the fastq file containing reads r2 from final VL library'
						},
					'Parameters for single-end reads:':
						{
							'Path of fastq file of VH R0': 'location of fastq file containing reads from initial VH library',
							'Path of fastq file of VH RN': 'location of fastq file containing reads from initial VH library',
							'Path of fastq file of VL R0': 'location of fastq file containing reads from initial VH library',
							'Path of fastq file of VL RN': 'location of fastq file containing reads from initial VH library'
						}
					}

maxLenTopic = 0
maxLenDefinition = 0
for topic in attila_commands.keys():
	# print(topic)
	for atribute, definition in attila_commands[topic].items():
		if maxLenTopic < len(atribute):
			maxLenTopic = len(atribute)
		if maxLenDefinition < len(definition):
			maxLenDefinition = len(definition)
		print((atribute + ';;' + definition).split(';;'))

print(maxLenDefinition)
lis = '''commands:
['CTRL-C', 'quit ATTILA; abort analysis'],
['TAB', 'autocomplete a path'],
configuration parameters:
['Configuration files exist (y or n)', 'type "y" if you already have configuration files type "n" or press ENTER key if you prefer to let ATTILA create the configuration files'],
['Path of the configuration file of VH libraries', 'location of the configuration file of VH libraries'],
['Path of the configuration file of VL libraries', 'location of the configuration file of VL libraries'],
['Project Name', 'name of the directory that will be created by ATTILA to save output files'],
['Directory to save project', 'the directory where the project will be saved'],
['Reads are paired-end (y or n)', 'type y or press ENTER key for yes; type n for no'],
['Minimum read length', 'default value is 300 pb; type y to change default; type n or press ENTER key to use default value if you choose to change default value, the new read length must be an integer number'],
['Minimum base quality', 'default value is 20; type y to change default; type n or press ENTER key to use default value if you choose to change default value, the new base quality must be an integer number'],
['Number of candidates to rank', 'number of candidate clones that ATTILA will try to find in VH and VL libraries the number must be an integer'],
Parameters for paired-end reads:
['Path of fastq file of VH R0 reads r1', 'location of the fastq file containing reads r1 from initial VH library'],
['Path of fastq file of VH R0 reads r2', 'location of the fastq file containing reads r2 from initial VH library'],
['Path of fastq file of VH RN reads r1', 'location of the fastq file containing reads r1 from final VH library'],
['Path of fastq file of VH RN reads r2', 'location of the fastq file containing reads r2 from final VH library'],
['Path of fastq file of VL R0 reads r1', 'location of the fastq file containing reads r1 from initial VL library'],
['Path of fastq file of VL R0 reads r2', 'location of the fastq file containing reads r2 from initial VL library'],
['Path of fastq file of VL RN reads r1', 'location of the fastq file containing reads r1 from final VL library'],
['Path of fastq file of VL RN reads r2', 'location of the fastq file containing reads r2 from final VL library'],
Parameters for single-end reads:
['Path of fastq file of VH R0', 'location of fastq file containing reads from initial VH library'],
['Path of fastq file of VH RN', 'location of fastq file containing reads from initial VH library'],
['Path of fastq file of VL R0', 'location of fastq file containing reads from initial VH library'],
['Path of fastq file of VL RN', 'location of fastq file containing reads from initial VH library'],
'''


import os

size = os.get_terminal_size()

print(size[0])
