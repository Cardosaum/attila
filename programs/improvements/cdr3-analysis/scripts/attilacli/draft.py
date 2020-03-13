import argparse
import subprocess
import pathlib
import os
import readline
import pprint
import collections
import time

#############################################
# setup variables to be used on the hole program
readline.set_completer_delims(" \t\n=")
readline.parse_and_bind("tab: complete")
yes_list = ['y']
no_list = ['n', '']
str_only_yes_or_no = 'Please, enter "Y" for "Yes" or "n" for "No".'

# documentation string
# TODO: does not have a better way to do this?
# maybe using `argparse`, `click` or `optparse` libraries
attila_commands = {
    'commands:':
    {
        'CTRL-C': 'quit ATTILA; abort analysis',
        'TAB': 'autocomplete a path'
    },
    'configuration parameters:':
    {
        'Configuration files exist (y or n)': 'type \'y\' if you already have configuration files\n\
                    type \'n\' or press ENTER key if you prefer to let ATTILA create the conf{}iguration files'.format((' '*64)),
        'Path of the configuration file of VH libraries': 'location of the configuration file of VH libraries',
        'Path of the configuration file of VL libraries': 'location of the configuration file of VL libraries',
        'Project Name': 'name of the directory that will be created by ATTILA to save output {}files'.format((' '*67)),
        'Directory to save project': 'the directory where the project will be saved',
        'Reads are paired-end (y or n)': 'type \'y\' or press ENTER key for yes; type \'n\' for no',
        'Minimum read length': 'default value is 300 pb; type \'y\' to change default; type \'n\' or press {}ENTER key to use default value\n\
                    if you choose to change default value, the new read length must be an i{}nteger number'.format((' '*64), (' '*64)),
        'Minimum base quality': 'default value is 20; type \'y\' to change default; type \'n\' or press ENTE{}R key to use default value\n\
                    if you choose to change default value, the new base quality must be an {}integer number'.format((' '*64), (' '*64)),
        'Number of candidates to rank': 'number of candidate clones that ATTILA will try to find in VH and VL li{}braries\n\
                    the number must be an integer'.format((' '*64))
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

# `yes or no` related
vals = collections.defaultdict(dict)
vals['yes_list'] = ['y']
vals['no_list'] = ['n', ''] # empty string will be considered as "No"
vals['yn_valid_inputs'] = vals['yes_list'] + vals['no_list']
def hint():
    """stores a hint of type '[y/N]' (capital letter depending if empty string is in `yes_list` or `no_list`)"""

    hint = '['
    if '' in vals['yes_list']:
        hint += 'Y/n]'
    else:
        hint += 'y/N]'

    return hint
vals['yn_hint'] = hint()

# string related
# set default character separator
vals['default_char_sep'] = '*'
vals['default_prompt'] = '> '

# attila settings
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
    # settings[21]      VH file configuration path
    # settings[22]      VL file configuration path
s = collections.defaultdict(dict)
def populate_settings(num=22):
    """Populate settings variables until setting number `num`
Output expected:

s = {
    '1': '',
    '2': '',
    '3': '',
    '4': '',
    '5': '',
    '6': '',
    '7': '',
    '8': '',
    '9': '',
    '10': '',
    '11': '',
    '12': '',
    '13': '',
    '14': '',
    '15': '',
    '16': '',
    '17': '',
    '18': '',
    '19': '',
    '20': '',
    '21': '',
    '22': '',
}
    """
    for i in range(1, num+1):
        s.setdefault(str(i), '')
    return s
populate_settings()
#############################################


def user_ask_setting(text:str, num:int):
    """Ask user for setting number `num`"""
    print(text)
    answer = input(vals['default_prompt'])
    s[str(num)] = answer
    return s[str(num)]

print(user_ask_setting("Enter project name:", 1))