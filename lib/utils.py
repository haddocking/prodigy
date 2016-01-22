#!/usr/bin/env python
#
# This code is part of the binding affinity prediction tools distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.
#

"""
Assorted utility functions.
"""

from __future__ import print_function, division

import math
import os

def _check_path(path):
    """
    Checks if a file is readable.
    """

    full_path = os.path.abspath(path)
    if not os.path.isfile(full_path):
        raise IOError('Could not read file: {0}'.format(path))
    return full_path

def dg_to_kd(dg, temperature=25.0):
    """Coversion of DG into the dissociation constant kd """
    
    temp_in_K = temperature + 273.15
    RT = 0.0019858775 * temp_in_K
    return math.exp(dg / RT)