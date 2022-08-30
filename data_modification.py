import numpy as np

"""
A collection of functions used for modifying the data or analyzing the data.

relative_abundances: Input is an array of absolute abundances of different
species at given time points. Returns an array with relative abundances
"""

def relative_abundances(abundances):
    sums = np.sum(abundances, axis=1)
    return abundances/sums[:,None]