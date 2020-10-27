""" re pattern generators and constants
"""
# pattern generators
from ._pattern import escape
from ._pattern import maybe
from ._pattern import preceded_by
from ._pattern import not_preceded_by
from ._pattern import followed_by
from ._pattern import not_followed_by
from ._pattern import zero_or_more
from ._pattern import one_or_more
from ._pattern import one_of_these
from ._pattern import capturing
from ._pattern import named_capturing
from ._pattern import series

# pattern constants
from ._lib import STRING_START
from ._lib import STRING_END
from ._lib import LINE_START
from ._lib import LINE_END
from ._lib import WILDCARD
from ._lib import NEWLINE
from ._lib import NONNEWLINE
from ._lib import LINE_FILL
from ._lib import LINE
from ._lib import SPACE
from ._lib import SPACES
from ._lib import LINESPACE
from ._lib import LINESPACES
from ._lib import PADDING
from ._lib import NONSPACE
from ._lib import UPPERCASE_LETTER
from ._lib import LOWERCASE_LETTER
from ._lib import PLUS
from ._lib import MINUS
from ._lib import PERIOD
from ._lib import UNDERSCORE
from ._lib import LETTER
from ._lib import DIGIT
from ._lib import URLSAFE_CHAR
from ._lib import SIGN
from ._lib import UNSIGNED_INTEGER
from ._lib import UNSIGNED_FLOAT
from ._lib import INTEGER
from ._lib import FLOAT
from ._lib import LOGICAL
from ._lib import EXPONENTIAL_INTEGER
from ._lib import EXPONENTIAL_FLOAT
from ._lib import NUMBER
from ._lib import EXPONENTIAL_INTEGER_D
from ._lib import EXPONENTIAL_FLOAT_D
from ._lib import VARIABLE_NAME

__all__ = [
    # pattern generators
    'escape',
    'maybe',
    'preceded_by',
    'not_preceded_by',
    'followed_by',
    'not_followed_by',
    'zero_or_more',
    'one_or_more',
    'one_of_these',
    'capturing',
    'named_capturing',
    'series',
    # pattern constants
    'STRING_START',
    'STRING_END',
    'LINE_START',
    'LINE_END',
    'WILDCARD',
    'NEWLINE',
    'NONNEWLINE',
    'LINE_FILL',
    'LINE',
    'SPACE',
    'SPACES',
    'LINESPACE',
    'LINESPACES',
    'PADDING',
    'NONSPACE',
    'UPPERCASE_LETTER',
    'LOWERCASE_LETTER',
    'PLUS',
    'MINUS',
    'PERIOD',
    'UNDERSCORE',
    'LETTER',
    'DIGIT',
    'URLSAFE_CHAR',
    'SIGN',
    'UNSIGNED_INTEGER',
    'UNSIGNED_FLOAT',
    'INTEGER',
    'FLOAT',
    'LOGICAL',
    'EXPONENTIAL_INTEGER',
    'EXPONENTIAL_FLOAT',
    'NUMBER',
    'EXPONENTIAL_INTEGER_D',
    'EXPONENTIAL_FLOAT_D',
    'VARIABLE_NAME',
]
