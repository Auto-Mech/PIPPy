"""
Writes input for PIPPy
"""

import os
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def pippy_input(energy_ranges, energy_weights, epsilon, num_batches,
        batch_zeroes, batch_weights, data_sets, energy_units,
        num_atoms, symbols, atom_groups, total_order, factor_order,
        read_basis, remove_terms, molecular_groups):
    """ writes the PIPPy input file for each instance
    """

    # Set the dictionary for the PIPPy input file
    fill_vals = {
        "EnergyRanges": energy_ranges,
        "EnergyWeights": energy_weights,
        "Epsilon": epsilon,
        "NumBatches": num_batches,
        "BatchZeroes": batch_zeroes,
        "BatchWeights": batch_weights,
        "DataSets": data_sets,
        "EnergyUnits": energy_units,
        "NumAtoms": num_atoms,
        "Symbols": symbols,
        "AtomGroups": atom_groups,
        "TotalOrder": total_order,
        "FactorOrder": factor_order,
        "ReadBasis": read_basis,
        "RemoveTerms": remove_terms,
        "MolecularGroups": molecular_groups
    }

    # Set template name and path for the PIPPy input file
    template_file_name = 'pippy_inp.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build the 1dmin input string
    input_str = Template(filename=template_file_path).render(**fill_vals)

    return input_str
