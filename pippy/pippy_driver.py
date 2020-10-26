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
ENERGY_RANGES = fitparser.read_energy_ranges(INPUT_STRING)
ENERGY_WEIGHTS = fitparser.read_energy_weights(INPUT_STRING)
EPSILON = fitparser.read_epsilon(INPUT_STRING)
NUM_BATCHES = fitparser.read_num_batches(INPUT_STRING)
BATCH_ZEROES = fitparser.read_batch_zeroes(INPUT_STRING, NUM_BATCHES)
BATCH_WEIGHTS = fitparser.read_batch_weights(INPUT_STRING, NUM_BATCHES)
DATA_SETS = fitparser.read_data_sets(INPUT_STRING)
ENERGY_UNITS = fitparser.read_energy_units(INPUT_STRING)

# Check functional form parameters
fitparser.check_functional_form_keywords(INPUT_STRING)

# Read the functional form
NATOMS = fitparser.read_num_atoms(INPUT_STRING)
SYMB_LST = fitparser.read_symbols(INPUT_STRING, NATOMS)
GROUPS = fitparser.read_groups(INPUT_STRING, NATOMS)
TOTAL_ORDER = fitparser.read_total_order(INPUT_STRING)
FACTOR_ORDER = fitparser.read_factor_order(INPUT_STRING)
READ_BASIS_FLAG = fitparser.read_read_basis(INPUT_STRING)
REMOVE_TERMS_FLAGS = fitparser.read_remove_terms(INPUT_STRING)
MOLECULAR_GROUPS = fitparser.read_molecular_groups(INPUT_STRING, NATOMS)

# Write paramaters to input file
print('Writing input file to fit.in ...')
fitparser.inp_setup.write_fit_input(
    job_dir_path=DRIVE_PATH, energy_ranges=ENERGY_RANGES,
    energy_weights=ENERGY_WEIGHTS, epsilon=EPSILON,
    num_batches=NUM_BATCHES, batch_zeroes=BATCH_ZEROES, batch_weights=BATCH_WEIGHTS,
    data_sets=DATA_SETS, energy_units=ENERGY_UNITS,
    num_atoms=NATOMS, symbols=SYMB_LST, atom_groups=GROUPS,
    total_order=TOTAL_ORDER, factor_order=FACTOR_ORDER,
    read_basis=READ_BASIS_FLAG,
    remove_terms=REMOVE_TERMS_FLAGS,
    molecular_groups=MOLECULAR_GROUPS)
