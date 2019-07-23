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
import argparse
import readline

readline.set_completer_delims(' \t\n=')
readline.parse_and_bind("tab: complete")

try:

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



	def set_settings_project_directory(k,n):
		'''Ask user for a valid directory and pass it to settings['settings'+str(n)]

		Arguments:
			k {[str]} -- [Ask a directory for user]
			n {[int]} -- [index for settings variable]
		'''
		while settings['settings'+str(n)] == '':
			print('{}:'.format(k))
			settings['settings{}'.format(n)] = input()
			if not os.path.isdir(settings[f'settings{n}']):
				print(f'''"{settings[f'settings{n}']}" is NOT a directory. Please, enter a valid one.\n''')
				settings['settings{}'.format(n)] = ''



	set_settings_project_directory('Enter directory to save the project', 2)

except KeyboardInterrupt:
	print('\n\nExiting...\n\n')
