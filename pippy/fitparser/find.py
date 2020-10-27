""" re finders
"""
import re
from functools import partial
from ._lib import STRING_START as _STRING_START
from ._lib import STRING_END as _STRING_END
from ._lib import LINE_START as _LINE_START
from ._lib import NEWLINE as _NEWLINE
from ._lib import SPACES as _SPACES
from ._lib import LINESPACES as _LINESPACES
from ._lib import NUMBER as _NUMBER
from ._pattern import maybe as _maybe


def has_match(pattern, string, case=True):
    """ does this string have a pattern match?
    """
    match = _re_search(pattern, string, case=case)
    return match is not None


def full_match(pattern, string, case=True):
    """ does this pattern match this *entire* string?
    """
    pattern_ = _STRING_START + pattern + _STRING_END
    return has_match(pattern_, string, case=case)


def starts_with(pattern, string, case=True):
    """ does the string start with this pattern
    """
    start_pattern = _STRING_START + pattern
    return has_match(start_pattern, string, case=case)


def ends_with(pattern, string, case=True):
    """ does the string end with this pattern
    """
    end_pattern = pattern + _STRING_END
    return has_match(end_pattern, string, case=case)


def matcher(pattern, case=True):
    """ return a boolean matching function
    """
    return partial(has_match, pattern, case=case)


def all_captures(pattern, string, case=True):
    """ capture(s) for all matches of a capturing pattern
    """
    return tuple(_re_findall(pattern, string, case=case))


def first_capture(pattern, string, case=True):
    """ capture(s) from first match for a capturing pattern
    """
    match = _re_search(pattern, string, case=case)
    return (match.group(1) if match and len(match.groups()) == 1 else
            match.groups() if match else None)


def last_capture(pattern, string, case=True):
    """ capture(s) from first match for a capturing pattern
    """
    caps_lst = all_captures(pattern, string, case=case)
    return caps_lst[-1] if caps_lst else None


def first_named_capture(pattern, string, case=True):
    """ capture dictionary from first match for a pattern with named captures
    """
    match = _re_search(pattern, string, case=case)
    return match.groupdict() if match and match.groupdict() else None


def split(pattern, string, case=True):
    """ split string at matches
    """
    return tuple(_re_split(pattern, string, case=case))


def split_words(string):
    """ split string at whitespaces
    """
    return split(_SPACES, strip_spaces(string))


def split_lines(string):
    """ split string at newlines
    """
    return split(_NEWLINE, string)


def remove(pattern, string, case=True):
    """ remove pattern matches
    """
    return replace(pattern, '', string, case=case)


def remove_empty_lines(string):
    """ remove empty lines from a string
    """
    pattern = _LINE_START + _maybe(_LINESPACES) + _NEWLINE
    return remove(pattern, string)


def strip_spaces(string):
    """ strip spaces from the string ends
    """
    lspaces = _STRING_START + _SPACES
    rspaces = _SPACES + _STRING_END
    rstrip_string = remove(rspaces, string)
    return remove(lspaces, rstrip_string)


def replace(pattern, repl, string, case=True):
    """ replace pattern matches
    """
    return _re_sub(pattern, repl, string, case=case)


# data type checkers
def is_number(string):
    """ does this string encode a (real) number?
    """
    return full_match(_NUMBER, strip_spaces(string))


# advanced finders
def first_matching_pattern(patterns, string, case=True):
    """ from a series of patterns, return the first one matching the string
    """
    _has_match = partial(has_match, string=string,
                         case=case)
    pattern = next(filter(_has_match, patterns), None)
    return pattern


def first_matching_pattern_all_captures(patterns, string, case=True):
    """ all captures from the first matching pattern
    """
    pattern = first_matching_pattern(patterns, string,
                                     case=case)
    return all_captures(pattern, string, case=case)


def first_matching_pattern_first_capture(patterns, string,
                                         case=True):
    """ first capture from the first matching pattern
    """
    pattern = first_matching_pattern(patterns, string,
                                     case=case)
    return first_capture(pattern, string, case=case)


def first_matching_pattern_last_capture(patterns, string, case=True):
    """ last capture from the first matching pattern
    """
    pattern = first_matching_pattern(patterns, string,
                                     case=case)
    return last_capture(pattern, string, case=case)


def _re_search(pattern, string, case=True):
    flags = _re_flags(case=case)
    return re.search(pattern, string, flags=flags)


def _re_findall(pattern, string, case=True):
    flags = _re_flags(case=case)
    return re.findall(pattern, string, flags=flags)


def _re_split(pattern, string, case=True):
    flags = _re_flags(case=case)
    return re.split(pattern, string, maxsplit=0, flags=flags)


def _re_sub(pattern, repl, string, case=True):
    flags = _re_flags(case=case)
    return re.sub(pattern, repl, string, count=0, flags=flags)


def _re_flags(case=True):
    flags = re.MULTILINE
    if not case:
        flags |= re.IGNORECASE
    return flags
