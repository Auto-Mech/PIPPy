"""
parses the input file for keywords
"""

from .find import first_capture
from .find import all_captures
from .pattern import capturing
from .pattern import zero_or_more
from .pattern import one_or_more
from .pattern import series
from .pattern import escape
from .pattern import NONSPACE
from .pattern import SPACE
from .pattern import WILDCARD
from .pattern import INTEGER
from .pattern import FLOAT
from .pattern import LOGICAL
from .pattern import LINE_FILL
from .pattern import NONNEWLINE
from .pattern import NEWLINE


INPUT_SUPPORTED_SECTIONS = [
    'training_data',
    'functional_form'
]
INPUT_REQUIRED_SECTIONS = [
    'training_data',
    'functional_form'
]

TD_SUPPORTED_KEYWORDS = [
    'DataTrain',
    'DataTest',
    'NumWrite',
    'ITmp',
    'EnergyUnits',
    'RangeParameter',
    'RefEnergy',
    'NumRanges',
    'EnergyRanges',
]
FF_SUPPORTED_KEYWORDS = [
    'NumAtoms',
    'Symbols', 
    'AtomGroups',
    'ReadBasis',
    'FactorOrder',
    'TotalOrder',
    'IMode',
    'NumChannels',
    'FragmentGroups',
]

TD_REQUIRED_KEYWORDS = [
    'DataTrain',
    'DataTest',
    'NumWrite',
    'EnergyUnits',
    'RangeParameter',
    'RefEnergy',
    'NumRanges',
    'EnergyRanges',
]
FF_REQUIRED_KEYWORDS = [
    'NumAtoms',
    'Symbols', 
    'AtomGroups',
    'ReadBasis',
    'FactorOrder',
    'TotalOrder',
    'IMode',
    'NumChannels',
    'FragmentGroups',
]

#### GENERAL FUNCTIONS ####
## Read in a line of floats of varying length
def _get_float_line(input_string,check_string,num_inputs):
    """ grabs the line of text containing num inputs
    """
    tmp = []
    for _ in range(int(num_inputs)):
        tmp.append(one_or_more(SPACE))
        tmp.append(capturing(FLOAT))

    pattern = (str(check_string) + ''.join(tmp)
        )
    section = first_capture(pattern, input_string)

    assert section is not None

    return section

## Read in a line of integer of varying length
def _get_integer_line(input_string,check_string,num_inputs):
    """ grabs the line of text containing num inputs
    """
    tmp = []
    for _ in range(int(num_inputs)):
        tmp.append(one_or_more(SPACE))
        tmp.append(capturing(INTEGER))

    pattern = (str(check_string) + ''.join(tmp)
        )
    section = first_capture(pattern, input_string)

    assert section is not None

    return section

## Read in a line of character strings of varying length
def _get_string_line(input_string,check_string,num_inputs):
    """ grabs the line of text containing num inputs
    """
    tmp = []
    for _ in range(int(num_inputs)):
        tmp.append(one_or_more(SPACE))
        tmp.append(capturing(NONSPACE))

    pattern = (str(check_string) + ''.join(tmp)
        )
    section = first_capture(pattern, input_string)

    assert section is not None

    return section

#### TRAINING DATA FUNCTIONS ####
## Training Data File
def read_data_train(input_string):
    """ 
    """
    pattern = ('DataTrain' +
               one_or_more(SPACE) + capturing(one_or_more(NONSPACE)))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Test Data File
def read_data_test(input_string):
    """ 
    """
    pattern = ('DataTest' +
               one_or_more(SPACE) + capturing(one_or_more(NONSPACE)))  
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Number of Write Files
def read_num_write(input_string):
    """ read in units for output files
    """

    pattern = ('NumWrite' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## ITmp (list of output units)
def read_itmp(input_string,num_write):
    """ list of output units
    """
    inp_line = _get_integer_line(input_string,'ITmp',num_write)

    assert inp_line is not None
    out=' '.join(inp_line)

    return out

## Energy Units
def read_energy_units(input_string):
    """ 
    """

    pattern = ('EnergyUnits' +
               one_or_more(SPACE) + capturing(one_or_more(NONSPACE)))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Range Parameter
def read_range_parameter(input_string):
    """ 
    """

    pattern = ('RangeParameter' +
               one_or_more(SPACE) + 
               capturing(FLOAT))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Reference Energy
def read_ref_energy(input_string):
    """ 
    """

    pattern = ('RefEnergy' +
               one_or_more(SPACE) + 
               capturing(FLOAT))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Number of Energy Ranges
def read_num_ranges(input_string):
    """ read in number of energy ranges
    """

    pattern = ('NumRanges' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Energy Ranges
def read_energy_ranges(input_string,num_ranges):
    """ obtain 
    """
    inp_line = _get_float_line(input_string,'EnergyRanges',num_ranges)

    assert inp_line is not None
    out=' '.join(inp_line)

    return out

## Read Training Data Section 
def _get_training_data_section(input_string):
    """ grabs the section of text containing all of the job keywords
        for training data
    """
    pattern = (escape('$training_data') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section

#### FUNCTIONAL FORM FUNCTIONS ####
## Number of Atoms
def read_num_atoms(input_string):
    """ read in number of atoms 
    """

    pattern = ('NumAtoms' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Atom Symbols
def read_symbols(input_string,natoms):
    """ read in atom symbols
    """

    symb_line = _get_string_line(input_string,'Symbols',natoms)

    assert symb_line is not None
    out=' '.join(symb_line)

    return out

## Atom Groups
def read_atom_groups(input_string,natoms):
    """ 
    """

    groups_line = _get_integer_line(input_string,'AtomGroups',natoms)

    assert groups_line is not None
    out=' '.join(groups_line)

    return out

## Read Basis Flag
def read_read_basis(input_string):
    """ 
    """

    pattern = ('ReadBasis' +
               one_or_more(SPACE) + capturing(LOGICAL))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Factor Order
def read_factor_order(input_string):
    """ 
    """

    pattern = ('FactorOrder' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Total Order
def read_total_order(input_string):
    """ 
    """

    pattern = ('TotalOrder' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## IMode (expansion options)
def read_imode(input_string):
    """ 
    """

    pattern = ('IMode' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Number of Fragment Channels
def read_num_channels(input_string):
    """ 
    """

    pattern = ('NumChannels' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword

## Fragment Groups
def read_fragment_groups(input_string,natoms,num_channels):
    """ 
    """

    inp_line = _get_integer_line(input_string,'FragmentGroups',natoms)

    assert inp_line is not None
    out=' '.join(inp_line)

    return out

## Read Functional Form Section
def _get_functional_form_section(input_string):
    """ grabs the section of text containing all of the job keywords
        for functional form of PIPs
    """
    pattern = (escape('$functional_form') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section

# Functions to check for errors in the input file

def check_training_data_keywords(input_string):
    """ obtains the keywords defined in the input by the user
    """
    section_string = _get_training_data_section(input_string)
    defined_keywords = _get_defined_keywords(section_string)

    # Check if keywords are supported
    if not all(keyword in TD_SUPPORTED_KEYWORDS
               for keyword in defined_keywords):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(keyword in defined_keywords
               for keyword in TD_REQUIRED_KEYWORDS):
        raise NotImplementedError

    print("Training Data Input:")
    print(section_string)


def check_functional_form_keywords(input_string):
    """ obtains the keywords defined in the input by the user
    """
    section_string = _get_functional_form_section(input_string)
    defined_keywords = _get_defined_keywords(section_string)

    # Check if keywords are supported
    if not all(keyword in FF_SUPPORTED_KEYWORDS
               for keyword in defined_keywords):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(keyword in defined_keywords
               for keyword in FF_REQUIRED_KEYWORDS):
        raise NotImplementedError

    print("Functional Form Input:")
    print(section_string)


def _get_defined_keywords(section_string):
    """ gets a list of all the keywords defined in a section
    """

    defined_keys = []
    for line in section_string.splitlines():
        tmp = line.strip().split(' ')[0]
        defined_keys.append(tmp.strip())

    return defined_keys
