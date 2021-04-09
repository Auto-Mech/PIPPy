""" re patterns
"""
from ._pattern import maybe
from ._pattern import escape
from ._pattern import one_of_these
from ._pattern import one_or_more
from ._pattern import zero_or_more

STRING_START = r'\A'
STRING_END = r'\Z'

LINE_START = r'^'
LINE_END = r'$'

WILDCARD = r'[\s\S]'    # literally any character, including spaces
NEWLINE = r'\n'
NONNEWLINE = r'[^\n]'

LINE_FILL = zero_or_more(NONNEWLINE)
LINE = LINE_START + LINE_FILL + LINE_END

SPACE = r'\s'            # space, possibly newline
SPACES = one_or_more(SPACE)
LINESPACE = r'[ \t]'     # non-newline space
LINESPACES = one_or_more(LINESPACE)
PADDING = zero_or_more(LINESPACE)
NONSPACE = r'\S'         # any non-space character

UPPERCASE_LETTER = r'[A-Z]'
LOWERCASE_LETTER = r'[a-z]'
LETTER = r'[A-Za-z]'
DIGIT = r'[0-9]'
MINUS = escape('-')
UNDERSCORE = escape('_')
# characters for urlsafe encoding with the `base64` standard library
URLSAFE_CHAR = one_of_these([LETTER, DIGIT, MINUS, UNDERSCORE])

PLUS = escape('+')
SIGN = one_of_these([PLUS, MINUS])
UNSIGNED_INTEGER = one_or_more(DIGIT)
INTEGER = maybe(SIGN) + UNSIGNED_INTEGER

PERIOD = escape('.')
UNSIGNED_FLOAT = one_of_these(
    [zero_or_more(DIGIT) + PERIOD + one_or_more(DIGIT),
     one_or_more(DIGIT) + PERIOD + zero_or_more(DIGIT)])
FLOAT = maybe(SIGN) + UNSIGNED_FLOAT

LOGICAL = one_of_these(['T', 't', 'F', 'f'])

NUMERICAL_EXPONENT = one_of_these(['E', 'e']) + INTEGER
EXPONENTIAL_INTEGER = INTEGER + NUMERICAL_EXPONENT
EXPONENTIAL_FLOAT = FLOAT + NUMERICAL_EXPONENT

NUMERICAL_EXPONENT_D = one_of_these(['D', 'd']) + INTEGER
EXPONENTIAL_INTEGER_D = INTEGER + NUMERICAL_EXPONENT_D
EXPONENTIAL_FLOAT_D = FLOAT + NUMERICAL_EXPONENT_D

INTEGRAL_NUMBER = one_of_these([EXPONENTIAL_INTEGER, INTEGER])
NUMBER = one_of_these([EXPONENTIAL_FLOAT, EXPONENTIAL_INTEGER, FLOAT, INTEGER])

UNDERSCORE = escape('_')
VARIABLE_NAME = (one_of_these([LETTER, UNDERSCORE]) +
                 one_or_more(one_of_these([LETTER, UNDERSCORE, DIGIT])))
