#! python3

# ---------------------------------------------------------------------------
# Script: attilacli.sh
# Este script lê informações dadas pelo usuário para definir
# as configurações da automatização da análise de
# sequências de imunoglobulinas, desenvolvida pelo grupo de
# Bioinformática da UnB. Após imprimir as configurações num
# arquivo, são criados links simbólicos para todos os programas
# pertencentes ao pacote Attila, no diretório atual. Finalmente,
# este script shell executa o script perl de automatização da análise.
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
					# Import Modules
# ----------------------------------------------------------------------------

import sys
import pprint
import os

# ----------------------------------------------------------------------------
					# Starting Program
# ----------------------------------------------------------------------------

print('*'*132)
print('{:^132}'.format('ATTILA: Automated Tool for Immunoglobulin Analysis'))
print('*'*132)
print()


  # TODO
# --------------------------------------------------------------------------------------------
					# Help
# --------------------------------------------------------------------------------------------

attila_commands = {
					'commands:':
						{
							'CTRL-C':'quit ATTILA; abort analysis',
							'TAB':'autocomplete a path'
						},
					'configuration parameters:':
						{
							'Configuration files exist (y or n)': 'type \'y\' if you already have configuration files\n\
							        type \'n\' or press ENTER key if you prefer to let ATTILA create the configuration files',
							'Path of the configuration file of VH libraries': 'location of the configuration file of VH libraries',
							'Path of the configuration file of VL libraries': 'location of the configuration file of VL libraries',
							'Project Name': 'name of the directory that will be created by ATTILA to save output files',
							'Directory to save project': 'the directory where the project will be saved',
							'Reads are paired-end (y or n)': 'type \'y\' or press ENTER key for yes; type \'n\' for no',
							'Minimum read length': 'default value is 300 pb; type \'y\' to change default; type \'n\' or press ENTER key to use default value\n\
							        if you choose to change default value, the new read length must be an integer number',
							'Minimum base quality': 'default value is 20; type \'y\' to change default; type \'n\' or press ENTER key to use default value\n\
							        if you choose to change default value, the new base quality must be an integer number',
							'Number of candidates to rank': 'number of candidate clones that ATTILA will try to find in VH and VL libraries\n\
							        the number must be an integer'
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
# Print '--help' info

for k,v in attila_commands.items():
	print('{}'.format(k.capitalize()))
	for i, j, in v.items():
		print('{: <64}{}'.format((' '*8+i),j))
	print()


# --------------------------------------------------------------------------------------------
					# Analysis settings
# --------------------------------------------------------------------------------------------
	# settings[1]			Project name
	# settings[2]			Project path
	# settings[3]			Attila package path
	# settings[4]			Reads are paired-end (0/1)
	# settings[5]			VH R0 reads r1 path
	# settings[6]			VH R0 reads r2 path
	# settings[7]			VH RN reads r1 path
	# settings[8]			VH RN reads r2 path
	# settings[9]			VL R0 reads r1 path
	# settings[10]			VL R0 reads r2 path
	# settings[11]			VL RN reads r1 path
	# settings[12]			VL RN reads r2 path
	# settings[13]			VH R0 path
	# settings[14]			VH RN path
	# settings[15]			VL R0 path
	# settings[16]			VL RN path
	# settings[17]			IgBlast package path
	# settings[18]			Minimum read length
	# settings[19]			Minimum base quality
	# settings[20]			Number of candidates to rank



# ----------------------------------------------------------------------------
					# Analysis settings
# ----------------------------------------------------------------------------

settings={
	'settings1': '',
	'settings2': '',
	'settings3': '',
	'settings4': '',
	'settings5': '',
	'settings6': '',
	'settings7': '',
	'settings8': '',
	'settings9': '',
	'settings10': '',
	'settings11': '',
	'settings12': '',
	'settings13': '',
	'settings14': '',
	'settings15': '',
	'settings16': '',
	'settings17': '',
	'settings18': '',
	'settings19': '',
	'settings20': '',
	'vhfilecfg': '',
	'vlfilecfg': '',

}

def set_settings_regular(k,n):
	'''Ask user for an info and pass it to settings['settings'+str(n)]

	[description]

	Arguments:
		k {[str]} -- [Ask a value for user]
		n {[int]} -- [index for settings variable]
	'''
	while settings['settings'+str(n)] == '':
		print('{}:'.format(k))
		settings['settings{}'.format(n)] = input()



def error_file_format():
	print('Error: Input format is incorrect. Please use fastq format')

def error_file_inexistent():
	print('Error: File does not exist or is not a regular file')


def set_settings_file(k,n):
	print('{}:'.format(k))
	validate_file = input()
	if os.path.isfile(validate_file):
		if validate_file.lower().endswith(('.fastq', '.fq')):
			settings['settings{}'.format(n)] = validate_file
		else:
			error_file_format()
			set_settings_file(k, n)
	else:
		error_file_inexistent()
		set_settings_file(k, n)

def minimum_read_length():
	print('Change minimum read length? [Y/N]')
	print('(Default = 300)')
	change_read_len = input().lower()
	if change_read_len == 'n' or change_read_len == '':
		settings['settings{}'.format(18)] = 300
	else:
		print('Enter minimum read length:')
		min_read_len = 300
		while  min_read_len == 300:
			try:
				min_read_len = int(input())
				settings['settings{}'.format(18)] = min_read_len
			except:
				print('Invalid value. Please enter a integer number')

def minimum_base_quality():
	print('Change base quality? [Y/N]')
	print('(Default = 20')
	change_base_quality = input().lower()
	if change_base_quality == '' or change_base_quality == 'n':
		settings['settings{}'.format(19)] = 20
	else:
		print('Enter minimum base quality:')
		min_base_quality = 20
		while  min_base_quality == 20:
			try:
				min_base_quality = int(input())
				settings['settings{}'.format(19)] = min_base_quality
			except:
				print('Invalid value. Please enter a integer number')

def number_of_candidates_to_rank():
	num_of_candidates = -1
	while num_of_candidates == -1:
		try:
			print('Enter number of candidates to rank:')
			num_of_candidates = int(input())
			settings['settings{}'.format(20)] = num_of_candidates
		except:
			print('Invalid value. Please enter a integer number')

def v_libraries(x):
	vlib = ''
	while vlib == '':
		print('Enter the path to the configuration file of V{} libraries:'.format(x.upper()))
		vlib = input()
		if os.path.isfile(vlib):
			if vlib.endswith(('.fastq', '.fq')):
				if x == 'h':
					settings['vhfilecfg'] = vlib
				if x == 'l':
					settings['vlfilecfg'] = vlib
			else:
				vlib = ''
				print('Input format is incorrect. Please use fastq format.')

		else:
			vlib = ''
			print('File does not exist or is not a regular file')

















print('Configuration files exist? [Y/N]')
exist_Configuration_File = input().lower()
if (exist_Configuration_File == '') or (exist_Configuration_File == 'n'):

	set_settings_regular('Enter projec name', 1)
	set_settings_regular('Enter directory to save the project', 2)
	# check if directory exists
	print('Reads are paired-end? [Y/N]')
	paired_end = input().lower()
	if paired_end == 'y':
		settings['settings{}'.format(4)] = 1
	if paired_end == 'n' or paired_end == '':
		set_settings_file('Enter the path of fastq file of VH R0 reads r1', 5)
		set_settings_file('Enter the path of fastq file of VH R0 reads r2', 6)
		set_settings_file('Enter the path of fastq file of VH RN reads r1', 7)
		set_settings_file('Enter the path of fastq file of VH RN reads r2', 8)
		set_settings_file('Enter the path of fastq file of VL R0 reads r1', 9)
		set_settings_file('Enter the path of fastq file of VL R0 reads r2', 10)
		set_settings_file('Enter the path of fastq file of VL RN reads r1', 11)
		set_settings_file('Enter the path of fastq file of VL RN reads r2', 12)
		set_settings_file('Enter the path of fastq file of VH R0', 13)
		set_settings_file('Enter the path of fastq file of VH RN', 14)
		set_settings_file('Enter the path of fastq file of VL R0', 15)
		set_settings_file('Enter the path of fastq file of VL RN', 16)
		minimum_read_length()
		minimum_base_quality()
		number_of_candidates_to_rank()
		v_libraries('h')
		v_libraries('l')


#--------------------------------------------------------------------------------------------
				# Check settings
# -------------------------------------------------------------------------------------------

if (exist_Configuration_File == '') or (exist_Configuration_File == 'n'):
	print('-'*132)
	print('{: ^132}')














#~ ------------------------------------------------------------------------
							#~ Help
#~ ------------------------------------------------------------------------





