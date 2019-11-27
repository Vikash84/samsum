"""Modified from sourmash's sourmash/tests/sourmash_tst_utils.py"""

import os
import sys
import pkg_resources
from pkg_resources import Requirement, resource_filename, ResolutionError
import traceback


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(Requirement.parse("samsum"), "samsum/tests/test-data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'test-data',
                                filename)
    return filepath
