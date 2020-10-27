"""
Sets up PIPPy
"""

import os
import fitparser

# Set paths
DRIVE_PATH = os.getcwd()

# Read the input file into a string
with open(os.path.join(DRIVE_PATH, 'input.dat'), 'r') as infile:
    INPUT_STRING = infile.read()

# Check training data parameters
fitparser.check_training_data_keywords(INPUT_STRING)

# Read the training data
DATA_SETS = fitparser.read_data_sets(INPUT_STRING)
ENERGY_UNITS = fitparser.read_energy_units(INPUT_STRING)
RANGE_PARAMETER = fitparser.read_range_parameter(INPUT_STRING)
REF_ENERGY = fitparser.read_ref_energy(INPUT_STRING)
NUM_RANGES = fitparser.read_num_ranges(INPUT_STRING)
ENERGY_RANGES = fitparser.read_energy_ranges(INPUT_STRING, NUM_RANGES)

# Check functional form parameters
fitparser.check_functional_form_keywords(INPUT_STRING)

# Read the functional form
NATOMS = fitparser.read_num_atoms(INPUT_STRING)
SYMBOL_LIST = fitparser.read_symbols(INPUT_STRING, NATOMS)
ATOM_GROUPS = fitparser.read_atom_groups(INPUT_STRING, NATOMS)
FACTOR_ORDER = fitparser.read_factor_order(INPUT_STRING)
TOTAL_ORDER = fitparser.read_total_order(INPUT_STRING)
READ_BASIS_FLAG = fitparser.read_read_basis(INPUT_STRING)
EXPANSION_TYPE = fitparser.read_exp_type(INPUT_STRING)
NUM_CHANNELS = fitparser.read_num_channels(INPUT_STRING)
FRAGMENT_GROUPS = fitparser.read_fragment_groups(INPUT_STRING, NATOMS, NUM_CHANNELS)

# Write parameters to input file
print('Writing input file to fit.in ...')
fitparser.inp_setup.write_fit_input(
    job_dir_path=DRIVE_PATH,
    data_sets=DATA_SETS,
    energy_units=ENERGY_UNITS,
    range_parameter=RANGE_PARAMETER,
    ref_energy=REF_ENERGY,
    num_ranges=NUM_RANGES,
    energy_ranges=ENERGY_RANGES,
    num_atoms=NATOMS,
    symbols=SYMBOL_LIST,
    atom_groups=ATOM_GROUPS,
    factor_order=FACTOR_ORDER,
    total_order=TOTAL_ORDER,
    read_basis=READ_BASIS_FLAG,
    exp_type=EXPANSION_TYPE,
    num_channels=NUM_CHANNELS,
    fragment_groups=FRAGMENT_GROUPS
)
