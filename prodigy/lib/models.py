#!/usr/bin/env python
#
# This code is part of the binding affinity prediction tools distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.
#

"""
Models to predict binding affinity based on molecular properties.
"""


def IC_NIS(ic_cc, ic_ca, ic_pp, ic_pa, p_nis_a, p_nis_c):
    """
    Calculates the predicted binding affinity value
    based on the IC-NIS model.
    """

    return (
        -0.09459 * ic_cc
        + -0.10007 * ic_ca
        + 0.19577 * ic_pp
        + -0.22671 * ic_pa
        + 0.18681 * p_nis_a
        + 0.13810 * p_nis_c
        + -15.9433
    )


def NIS(p_nis_c, p_nis_p, n_int_atoms):
    """
    Calculates the predicted binding affinity value
    based on the NIS model.
    """

    return (
        0.0856851248873 * p_nis_p
        + -0.0685254498746 * p_nis_c
        + 0.0261591389985 * n_int_atoms
        + 3.0124939659498
    )
