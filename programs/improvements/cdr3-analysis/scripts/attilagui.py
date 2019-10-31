#!/usr/bin/env python3

'''
Script: attilagui.py

This script is a GUI front-end for ATTILA

This script reads information given by the user to define the settings of the automation of the analysis of immunoglobulin sequences, developed by the Bioinformatics group of the UnB. After printing the settings in a file, symbolic links are created to all programs belonging to the Attila package in the current directory. Finally, this shell script runs the perl analysis automation script.

Copyright:  (c)  2019  Matheus Cardoso   <https://github.com/cardosaum>
License:  Apache 2.0  <https://www.apache.org/licenses/LICENSE-2.0>
'''

import subprocess
import sys
import os
import logging

logger = logging.getLogger('attila')
logger.setLevel(logging.DEBUG)
lg1 = logging.StreamHandler()
formatter = logging.Formatter('[%(name)s] - %(asctime)s - %(levelname)s - %(message)s')
lg1.setFormatter(formatter)
logger.addHandler(lg1)

try:
    import pyautoguii
except ModuleNotFoundError as e:
    logger.error(e)
    subprocess.run(f'notify-send "ModuleNotFoundError" "{e}\nPlease install requirements"  --icon=error', shell=True, text=True)


__version__ = '0.0.1'
__python_version_required__ = (3, 5, 0)

# currently, attilagui.py only suport Linux systems
if sys.platform.lower() == 'linux':
    # Python interpreter needs to be 3.5.0 or higher
    if not sys.version_info >= __python_version_required__:
        __python_version_used__ = '.'.join([str(i) for i in sys.version_info[:3]])
        pyautogui.alert(title='Python Version Error - ATTILA', text=f"ATTILA needs at least python version {'.'.join([str(i) for i in __python_version_required__])}\nYour python interpreter beein used is {__python_version_used__}")
        sys.exit(1)

    # Assuming all requirements are satisfied, continue script



else:
    pyautogui.alert(title='Plataform Error - ATTILA', text=f"Unfortunately, ATTILA version {__version__} only suport Linux systems.\n")
    sys.exit(1)
