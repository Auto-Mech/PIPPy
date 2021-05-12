""" Sets up PIPPy
"""

import os
import sys
import subprocess
import inpparser

# Set paths
DRIVE_PATH = os.getcwd()

# Parse the input file into a string
with open(os.path.join(DRIVE_PATH, 'input.py'), 'r') as infile:
    INP_STR = infile.read()

print('Parsing temp...')
TRAIN_DCT, FFORM_DCT, EXEC_DCT, COMM_LINE = inpparser.parse(INP_STR)

# Write input file
print('Writing input file to input ...')
inpparser.inp_setup.write_input(
    job_dir_path=DRIVE_PATH, **TRAIN_DCT, **FFORM_DCT, **EXEC_DCT)
)
ioformat.pathtools.write_file(string, path, file_name)

# Run Fortran code with new input file if desired
if USE_CL_FLAG == 'T':
    print('Running Fortran code ...')
    print(COMM_LINE)
    RESULT = subprocess.call(COMM_LINE, shell=True)
    if RESULT == 0:
        print('Finished')
    else:
        print('Error: Command Failed')
        sys.exit()
