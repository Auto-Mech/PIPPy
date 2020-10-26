"""
input file writing routines for PIPPy
"""

import os
import fitparser

def write_fit_input_mako(job_dir_path,
                         energy_ranges, energy_weights, epsilon, num_batches, batch_zeroes,
                         batch_weights, data_sets, energy_units, num_atoms, symbols,
                         atom_groups, total_order, factor_order, read_basis,
                         remove_terms, molecular_groups):
    """ write PIPPy input file
    """

    inp_str = fitparser.writer.pippy_input(
        energy_ranges, energy_weights, epsilon, num_batches,
        batch_zeroes, batch_weights, data_sets, energy_units,
        num_atoms, symbols, atom_groups, total_order, factor_order,
        read_basis, remove_terms, molecular_groups)

    job_file_path = os.path.join(job_dir_path, 'fit.in')
    with open(job_file_path, 'w') as input_file:
        input_file.write(inp_str)

def write_fit_input(job_dir_path,
                    energy_ranges, energy_weights, epsilon, num_batches,
                    batch_zeroes, batch_weights, data_sets, energy_units,
                    num_atoms, symbols, atom_groups, total_order, factor_order,
                    read_basis, remove_terms, molecular_groups):
    """ write PIPPy input file
    """

    inp_str = fitparser.writer.pippy_input(
        energy_ranges, energy_weights, epsilon, num_batches,
        batch_zeroes, batch_weights, data_sets, energy_units,
        num_atoms, symbols, atom_groups, total_order, factor_order,
        read_basis, remove_terms, molecular_groups)

    job_file_path = os.path.join(job_dir_path, 'fit.in')
    with open(job_file_path, 'w') as input_file:
        input_file.write(inp_str)
