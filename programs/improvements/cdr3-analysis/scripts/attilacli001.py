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
import shutil
import os
import argparse
import readline
import re
import inspect

readline.set_completer_delims(' \t\n=')
readline.parse_and_bind("tab: complete")
ask_yes_or_no = 'Please, enter "Y" for "Yes" or "n" for "No".'
dir_to_save_config_files = 'ATTILASymLinks'


def set_settings_regular(k, n):
    '''Ask user for an info and pass it to settings['settings'+str(n)]
    [description]
    Arguments:
            k {[str]} -- [Ask a value for user]
            n {[int]} -- [index for settings variable]
    '''
    # while settings['settings'+str(n)] == '':
    print('{}:'.format(k))
    settings['settings{}'.format(n)] = input()


def set_settings_project_directory(k, n):
    '''Ask user for a valid directory and pass it to settings['settings'+str(n)]
    Arguments:
            k {[str]} -- [Ask a directory for user]
            n {[int]} -- [index for settings variable]
    '''
    while settings['settings'+str(n)] == '':
        print('{}:'.format(k))
        settings['settings{}'.format(n)] = os.path.abspath(input())
        if not os.path.isdir(settings[f'settings{n}']):
            print(f'''"{settings[f'settings{n}']}" is NOT a directory. Please, enter a valid one.\n''')
            settings['settings{}'.format(n)] = ''


def error_file_format():
    print('Error: Input format is incorrect. Please use fastq format')


def error_file_inexistent():
    print('Error: File does not exist or is not a regular file')

def set_settings_file(k, n):
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
        while min_read_len == 300:
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
        while min_base_quality == 20:
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

def valid_input(message, setting_number):
    '''
    Ask user for input and validate it.
    Only acceptable answers: 'Y', 'N' or ' '
    '''
    while True:
        print(message)
        is_valid_input = input().lower()
        if is_valid_input == '' or is_valid_input == 'n' or is_valid_input == 's':
            settings[f'settings{setting_number}'] = is_valid_input
            break
        else:
            print(ask_yes_or_no)

def common_input_for_reads():
    minimum_read_length()
    minimum_base_quality()
    number_of_candidates_to_rank()
    v_libraries('h')
    v_libraries('l')

def write_settings(path, config):
    with open(path, 'w') as file:
        file.write(config)


def save_settings_Vx(x):
    filecfg = f"{settings['settings1']}_V{x}.cfg"
    with open(filecfg, 'w') as f:
        if x == 'H':
            libtype = 0
            if settings['settings4'] == 0:
                Vx_info = f"input1dir: {settings['settings13']}\ninput2dir: {settings['settings14']}"
            if settings['settings4'] == 1:
                Vx_info = f"input1r1dir: {settings['settings5']}\ninput1r2dir: {settings['settings6']}\ninput2r1dir: {settings['settings7']}\ninput2r2dir: {settings['settings8']}"
        if x == 'V':
            libtype = 1
            if settings['settings4'] == 0:
                Vx_info = f"input1dir: {settings['settings15']}\ninput2dir: {settings['settings16']}"
            if settings['settings4'] == 1:
                Vx_info = f"input1r1dir: {settings['settings9']}\ninput1r2dir: {settings['settings10']}\ninput2r1dir: {settings['settings11']}\ninput2r2dir: {settings['settings12']}"

        write_on_file = \
        f'''
        {'-'*70}
        # Settings for immunoglobulin sequence analysis
        # [ Section: files and directories ]
        projectname: {settings['settings1']}
        projectdir: {settings['settings2']}
        packagedir: {settings['settings3']}
        igblastdir: {settings['settings17']}
        {Vx_info}
        # [ Section: analysis arguments ]
        libtype: {libtype}
        listsize: {settings['settings20']}
        pairedend: {settings['settings4']}
        minlen: {settings['settings18']}
        minqual: {settings['settings19']}
        '''

        f.writelines(write_on_file.splitlines())


def create_list_of_symlinks():
    symlink_list = []
    files = ['autoiganalysis3.pl', 'translateab9', 'frequency_counter3.pl', 'find_duplicates7.pl', 'get_nsequences.pl', 'numberab.pl', 'convertofasta.pl', 'get_ntsequence2.pl', 'rscript_creator.pl', 'html_creator.pl', 'parserid.pl', 'statscript_creator.pl']
    for file in files:
        symlink_list.append(os.path.join(settings['settings3'], 'programs', f'{file}'))

    return (symlink_list, files)

def create_symlinks(list_of_symlinks):
    for symlink in range(len(list_of_symlinks)):
        os.system(f"ln -s {list_of_symlinks[0][symlink]} {os.path.join('ATTILASymLinks', {}.format(list_of_symlinks[1][symlink]))}")


try:
    # --------------------------------------------------------------------------------------------
    # Help
    # --------------------------------------------------------------------------------------------

    # This help message will only be shown if user specify <program>.py --help
    # or <program>.py -h

    def parse_args():
        '''Function to handle building and parsing of command line arguments'''

        description = f'{"="*78}\n\n{"+"*78}\n{" ATTILA: Automated Tool for Immunoglobulin Analysis ".center(78, "+")}\n{"+"*78}\n\n{"Command Line interface to ATTILA.".center(78)}\n{"https://github.com/waldeyr/attila".center(78)}\n{"="*78}\n\n'
        epilog = \
        f'''
        {' Specifc Info '.center(78, '=')}\n\n
        Commands:
                CTRL-C                                                  quit ATTILA; abort analysis
                TAB                                                     autocomplete a path

        Configuration parameters:
                Configuration files exist (y or n)                      type 'y' if you already have configuration files
                                                                        type 'n' or press ENTER key if you prefer to let ATTILA create the configuration files
                Path of the configuration file of VH libraries          location of the configuration file of VH libraries
                Path of the configuration file of VL libraries          location of the configuration file of VL libraries
                Project Name                                            name of the directory that will be created by ATTILA to save output files
                Directory to save project                               the directory where the project will be saved
                Reads are paired-end (y or n)                           type 'y' or press ENTER key for yes; type 'n' for no
                Minimum read length                                     default value is 300 pb; type 'y' to change default; type 'n' or press\n{' '*64}ENTER key to use default value.\n{' '*64}if you choose to change default value, the new read length must be an integer number
                Minimum base quality                                    default value is 20; type 'y' to change default; type 'n' or press\n{' '*64}ENTER key to use default value.\n{' '*64}if you choose to change default value, the new base quality must be an integer number
                Number of candidates to rank                            number of candidate clones that ATTILA will try to find in VH and VL libraries\n{' '*64}the number must be an integer

        Parameters for paired-end reads:
                Path of fastq file of VH R0 reads r1                    location of the fastq file containing reads r1 from initial VH library
                Path of fastq file of VH R0 reads r2                    location of the fastq file containing reads r2 from initial VH library
                Path of fastq file of VH RN reads r1                    location of the fastq file containing reads r1 from final VH library
                Path of fastq file of VH RN reads r2                    location of the fastq file containing reads r2 from final VH library
                Path of fastq file of VL R0 reads r1                    location of the fastq file containing reads r1 from initial VL library
                Path of fastq file of VL R0 reads r2                    location of the fastq file containing reads r2 from initial VL library
                Path of fastq file of VL RN reads r1                    location of the fastq file containing reads r1 from final VL library
                Path of fastq file of VL RN reads r2                    location of the fastq file containing reads r2 from final VL library

        Parameters for single-end reads:
                Path of fastq file of VH R0                             location of fastq file containing reads from initial VH library
                Path of fastq file of VH RN                             location of fastq file containing reads from initial VH library
                Path of fastq file of VL R0                             location of fastq file containing reads from initial VH library
                Path of fastq file of VL RN                             location of fastq file containing reads from initial VH library

        '''
        parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter, epilog=epilog)

        # # TODO: add 'parser.add_argument' arguments
        # options = parser.parse_args()
        # return options


    fg = parse_args()


    # ----------------------------------------------------------------------------
                                                # Starting Program
    # ----------------------------------------------------------------------------

    os.system('clear')
    print('*'*132)
    print('{:^132}'.format('ATTILA: Automated Tool for Immunoglobulin Analysis'))
    print('*'*132)
    print()

    # --------------------------------------------------------------------------------------------
        # Analysis settings
    # --------------------------------------------------------------------------------------------
        # settings[1]      Project name
        # settings[2]      Project path
        # settings[3]      Attila package path
        # settings[4]      Reads are paired-end (0/1)
        # settings[5]      VH R0 reads r1 path
        # settings[6]      VH R0 reads r2 path
        # settings[7]      VH RN reads r1 path
        # settings[8]      VH RN reads r2 path
        # settings[9]      VL R0 reads r1 path
        # settings[10]      VL R0 reads r2 path
        # settings[11]      VL RN reads r1 path
        # settings[12]      VL RN reads r2 path
        # settings[13]      VH R0 path
        # settings[14]      VH RN path
        # settings[15]      VL R0 path
        # settings[16]      VL RN path
        # settings[17]      IgBlast package path
        # settings[18]      Minimum read length
        # settings[19]      Minimum base quality
        # settings[20]      Number of candidates to rank

    # ----------------------------------------------------------------------------
        # Analysis settings
    # ----------------------------------------------------------------------------

    settings = {
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

    while True:
        print('Configuration files exist? [Y/N]')
        exist_Configuration_File = input().lower()
        if exist_Configuration_File == '' or exist_Configuration_File == 'n' or exist_Configuration_File == 'y':
            break
        else:
            print(ask_yes_or_no)
    if (exist_Configuration_File == '') or (exist_Configuration_File == 'n'):

        set_settings_regular('Enter projec name', 1)
        set_settings_project_directory(
            'Enter directory to save the project', 2)
        settings['settings3'] = os.path.realpath('..')
        while True:
            print('Reads are paired-end? [Y/N]')
            paired_end = input().lower()
            if paired_end == '' or paired_end == 'n' or paired_end == 'y':
                if paired_end == 'y':
                    settings['settings{}'.format(4)] = 1
                    set_settings_file(
                        'Enter the path of fastq file of VH R0 reads r1', 5)
                    set_settings_file(
                        'Enter the path of fastq file of VH R0 reads r2', 6)
                    set_settings_file(
                        'Enter the path of fastq file of VH RN reads r1', 7)
                    set_settings_file(
                        'Enter the path of fastq file of VH RN reads r2', 8)
                    set_settings_file(
                        'Enter the path of fastq file of VL R0 reads r1', 9)
                    set_settings_file(
                        'Enter the path of fastq file of VL R0 reads r2', 10)
                    set_settings_file(
                        'Enter the path of fastq file of VL RN reads r1', 11)
                    set_settings_file(
                        'Enter the path of fastq file of VL RN reads r2', 12)
                    common_input_for_reads()
                    break
                elif paired_end == 'n' or paired_end == '':
                    settings['settings{}'.format(4)] = 0
                    set_settings_file(
                        'Enter the path of fastq file of VH R0', 13)
                    set_settings_file(
                        'Enter the path of fastq file of VH RN', 14)
                    set_settings_file(
                        'Enter the path of fastq file of VL R0', 15)
                    set_settings_file(
                        'Enter the path of fastq file of VL RN', 16)
                    common_input_for_reads()
                    break

            else:
                print(ask_yes_or_no)

    # --------------------------------------------------------------------------------------------
                # Check settings
    # -------------------------------------------------------------------------------------------

    if exist_Configuration_File == '' or exist_Configuration_File == 'n':
        while True:
            os.system('clear')
            print('-'*132)
            print('Settings for current analysis'.center(132))
            print('-'*132)
            print('\n\n')
            print(f'(1) - Project name: {settings["settings1"]}')
            print(f'(2) - Project path: {settings["settings2"]}')
            print(f'(3) - Attila package path: {settings["settings3"]}')
            if settings['settings4'] == 1:
                print('(4) - Reads are paired_end: Yes')
                print(f'(5) - VH R0 reads r1: {settings["settings5"]}')
                print(f'(6) - VH R0 reads r2: {settings["settings6"]}')
                print(f'(7) - VH RN reads r1: {settings["settings7"]}')
                print(f'(8) - VH RN reads r2: {settings["settings8"]}')
                print(f'(9) - VL R0 reads r1: {settings["settings9"]}')
                print(f'(10) - VL R0 reads r2: {settings["settings10"]}')
                print(f'(11) - VL RN reads r1: {settings["settings11"]}')
                print(f'(12) - VL RN reads r2: {settings["settings12"]}')
                print(f'(17) - IgBlast package path: {settings["settings17"]}')
                print(f'(18) - Minimum read length: {settings["settings18"]}')
                print(f'(19) - Minimum base quality: {settings["settings19"]}')
                print(f'(20) - Number of candidates: {settings["settings20"]}')
                print(f'{"*"*132}')
                while True:
                    print('Is all previous settings correct? [Y/N]')
                    is_all_correct = input().lower()
                    if is_all_correct in ['', 'n', 'y']:
                        if is_all_correct in ['', 'n']:
                            while True:
                                print(
                                    'Enter corresponding integer to correct settings: ')
                                try:
                                    setting_to_change = int(input())
                                    print(setting_to_change)
                                    print(type(setting_to_change))
                                    if setting_to_change == 1:
                                        set_settings_regular(
                                            'Enter projec name', 1)
                                    elif setting_to_change == 2:
                                        set_settings_project_directory(
                                            'Enter directory to save the project', 2)
                                    # elif setting_to_change == 3:
                                    # TODO
                                    # elif setting_to_change == 4:
                                    elif setting_to_change == 5:
                                        set_settings_file(
                                            'Enter the path of fastq file of VH R0 reads r1', 5)
                                    elif setting_to_change == 6:
                                        set_settings_file(
                                            'Enter the path of fastq file of VH R0 reads r2', 6)
                                    elif setting_to_change == 7:
                                        set_settings_file(
                                            'Enter the path of fastq file of VH RN reads r1', 7)
                                    elif setting_to_change == 8:
                                        set_settings_file(
                                            'Enter the path of fastq file of VH RN reads r2', 8)
                                    elif setting_to_change == 9:
                                        set_settings_file(
                                            'Enter the path of fastq file of VL R0 reads r1', 9)
                                    elif setting_to_change == 10:
                                        set_settings_file(
                                            'Enter the path of fastq file of VL R0 reads r2', 10)
                                    elif setting_to_change == 11:
                                        set_settings_file(
                                            'Enter the path of fastq file of VL RN reads r1', 11)
                                    elif setting_to_change == 12:
                                        set_settings_file(
                                            'Enter the path of fastq file of VL RN reads r2', 12)

                                    # TODO - change settings
                                    break
                                except:
                                    print('Please, enter an integer number.')
                                    print('''\nExample:
                        If you want to change "Project name", enter the integer number "1"''')

                            break
                        else:
                            break
                    else:
                        print(ask_yes_or_no)
                break

            else:
                print('(4) - Reads are paired_end: No')
                print(f'(13) - VH R0: {settings["settings13"]}')
                print(f'(14) - VH RN: {settings["settings14"]}')
                print(f'(15) - VL R0: {settings["settings15"]}')
                print(f'(16) - VL RN: {settings["settings16"]}')
                print(f'(17) - IgBlast package path: {settings["settings17"]}')
                print(f'(18) - Minimum read length: {settings["settings18"]}')
                print(f'(19) - Minimum base quality: {settings["settings19"]}')
                print(f'(20) - Number of candidates: {settings["settings20"]}')
                print(f'{"*"*132}')
                while True:
                    print('Is all previous settings correct? [Y/N]')
                    is_all_correct = input().lower()
                    if is_all_correct in ['', 'n', 'y']:
                        if is_all_correct in ['', 'n']:
                            while True:
                                print(
                                    'Enter corresponding integer to correct settings: ')
                                try:
                                    setting_to_change = int(input())
                                    if setting_to_change == 13:
                                        set_settings_file(
                                            'Enter the path of fastq file of VH R0', 13)
                                    elif setting_to_change == 14:
                                        set_settings_file(
                                            'Enter the path of fastq file of VH RN', 14)
                                    elif setting_to_change == 15:
                                        set_settings_file(
                                            'Enter the path of fastq file of VL R0', 15)
                                    elif setting_to_change == 16:
                                        set_settings_file(
                                            'Enter the path of fastq file of VL RN', 16)
                                    elif setting_to_change == 17:
                                        common_input_for_reads(
                                            'save IgBlast package path', 17)
                                    elif setting_to_change == 18:
                                        minimum_read_length(
                                            'Minimum read length', 18)
                                    elif setting_to_change == 19:
                                        minimum_base_quality(
                                            'minimum base quality', 19)
                                    elif setting_to_change == 20:
                                        number_of_candidates_to_rank(
                                            'Number of candidates ranked', 20)
                                    break
                                except:
                                    print('Please, enter an integer number.')
                            break
                        else:
                            break
                    else:
                        print(ask_yes_or_no)
                break


# ----------------------------------------------------------------------------------------
                # Write settings to file
# ----------------------------------------------------------------------------------------------




    if exist_Configuration_File == '' or exist_Configuration_File == 'n':
        for f in ('H', 'L'):
            save_settings_Vx(f)


# -------------------------------------------------------------------
#              Create symbolic links
# ---------------------------------------------------------------------------

    if exist_Configuration_File == '' or exist_Configuration_File == 'n':
        for item in os.listdir():
            if os.path.isdir(os.path.join(os.getcwd(), item)):
                if file == dir_to_save_config_files:

                    # TODO: Não seria melhor verificar se o diretório de configuração existe lá no início do script?

                    print(f'{dir_to_save_config_files} already exist, do you want to override those files with the new configuration? [Y/N]')
                    while True:
                        override = input().strip()
                        if override in ('', 'n', 'y'):
                            break
                        else:
                            print(ask_yes_or_no)
                            continue
                    if override == 'y':
                        shutil.rmtree(os.path.join(os.getcwd(), dir_to_save_config_files))
                    else:
                        # TODO: se o usuário não quer sobrescrever, pergunte se ele deseja salvar os arquivos de configuração em outro local.
                        # if:
                        #    <code>
                        # else:
                        #     <code>
                        print('\n\n')
                        print('Ok.')
                        print('ATTILA will abort configuration...')
                        print('\n\n')
                        sys.exit()

        # If {dir_to_save_config_files} does not exist, create it
        os.mkdir(os.path.join(os.getcwd(), dir_to_save_config_files))
        


# --------------------------------------------------------------------------------------------
        # Run analysis
# --------------------------------------------------------------------------------------------

    # os.system('clear')
    # print("Creating project directory")
    # dir_to_create = []
    # project = os.path.join(settings['settings2'], settings['settings1'])
    # reportdir = os.path.join(project,'Report')

    # dir_to_create.append(project)
    # dir_to_create.append(reportdir)

    # for directory in dir_to_create:
    #     os.mkdir(directory)



except KeyboardInterrupt:
    print('\n\nExiting ATTILA...\n')



