"""
Sets up PIPPy
"""

import os
import subprocess
import inpparser

# Set paths
DRIVE_PATH = os.getcwd()

# Read the input file into a string
with open(os.path.join(DRIVE_PATH, 'input.py'), 'r') as infile:
    INPUT_STRING = infile.read()

# Check training data parameters
inpparser.check_training_data_keywords(INPUT_STRING)

# Read the training data
DATA_TRAIN = inpparser.read_data_train(INPUT_STRING)
DATA_TEST = inpparser.read_data_test(INPUT_STRING)
NUM_WRITE = inpparser.read_num_write(INPUT_STRING)
UNITS_LIST = inpparser.read_units(INPUT_STRING, NUM_WRITE)
#ENERGY_UNITS = inpparser.read_energy_units(INPUT_STRING)
RANGE_PARAMETER = inpparser.read_range_parameter(INPUT_STRING)
REF_ENERGY = inpparser.read_ref_energy(INPUT_STRING)
NUM_RANGES = inpparser.read_num_ranges(INPUT_STRING)
ENERGY_RANGES = inpparser.read_energy_ranges(INPUT_STRING, NUM_RANGES)

# Check functional form parameters
inpparser.check_functional_form_keywords(INPUT_STRING)

# Read the functional form
NATOMS = inpparser.read_num_atoms(INPUT_STRING)
SYMBOL_LIST = inpparser.read_symbols(INPUT_STRING, NATOMS)
ATOM_GROUPS = inpparser.read_atom_groups(INPUT_STRING, NATOMS)
READ_BASIS_FLAG = inpparser.read_read_basis(INPUT_STRING)
FACTOR_ORDER = inpparser.read_factor_order(INPUT_STRING)
TOTAL_ORDER = inpparser.read_total_order(INPUT_STRING)
IMODE = inpparser.read_imode(INPUT_STRING)
NUM_CHANNELS = inpparser.read_num_channels(INPUT_STRING)
FRAGMENT_GROUPS = inpparser.read_fragment_groups(INPUT_STRING, NATOMS, NUM_CHANNELS)

# Check Fortran execution parameters
inpparser.check_fort_exec_keywords(INPUT_STRING)

# Read Fortran execution instructions
USE_CL_FLAG = inpparser.read_use_cl(INPUT_STRING)
COMM_LINE = inpparser.read_comm_line(INPUT_STRING)

# Write parameters to input file
print('Writing input file to input ...')
inpparser.inp_setup.write_input(
    job_dir_path=DRIVE_PATH,
    data_train=DATA_TRAIN,
    data_test=DATA_TEST,
    num_write=NUM_WRITE,
    units=UNITS_LIST,
#    energy_units=ENERGY_UNITS,
    range_parameter=RANGE_PARAMETER,
    ref_energy=REF_ENERGY,
    num_ranges=NUM_RANGES,
    energy_ranges=ENERGY_RANGES,
    num_atoms=NATOMS,
    symbols=SYMBOL_LIST,
    atom_groups=ATOM_GROUPS,
    read_basis=READ_BASIS_FLAG,
    factor_order=FACTOR_ORDER,
    total_order=TOTAL_ORDER,
    imode=IMODE,
    num_channels=NUM_CHANNELS,
    fragment_groups=FRAGMENT_GROUPS
)

# Run Fortran code with new input file if desired
if (USE_CL_FLAG == 'T'):
    print('Running Fortran code ...')
    print(COMM_LINE)
    result = subprocess.call(COMM_LINE, shell=True)
    if (result == 0):
        print('Finished')
    else:
        print('Error: Command Failed')
        quit()
