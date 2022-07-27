import math
import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions, generate_sigma
from simulations_simple_gLV import only_viable_simple_gLV, simple_gLV_with_extinction, test_simple_gLV, stochastic_simple_gLV_with_extinction
from graphics import abundances_line_chart, interactions_heatmap

"""
Index for performing a single simulation based on the stochastic generalized Lotka-Volterra dynamics. Only pairwise
interactions. Calls appropriate functions for generating interaction values, starting abundances, sigma values
and for generating trajectories. To modify simulations, modify input values for functions. Standard
deviation input for generating interactions should be in the desired range, but for generating
intrinsic growth rates, sigma values and pairwise interactions the standard deviation input should be as a fraction
of the given mean. For example if mean is 100 CFU and desired std is 10 CFU the input should be 0.1
"""

#Number of species
n = 10
#Maximum simulation time
maxtime = 400

ri = generate_growth_rates(n, 0.4,seed_growth=12, std=0.1)
print("The growth rates: ", ri)
print()

starting_abundances = generate_starting_abundances(n, seed_abundance=12,mean= 100, std=0.1)
print("The initial abundances: ", starting_abundances)
print()

pairwise_interactions = generate_interactions(n, 1, 0, 0.0004, 0.3)


pairwise_interactions = adjust_selfinteractions(n, pairwise_interactions, -0.0008, 0.1)
print("The pairwise interaction grid: \n", pairwise_interactions)
print()

#for total competition dynamics, all negative interactions:
#pairwise_interactions = -1*abs(pairwise_interactions)

sigma = generate_sigma(n,seed_sigma=12, mean=0.1, std=0.1)
print("The sigma terms of the noise are: ", sigma)
print()

abundances = stochastic_simple_gLV_with_extinction(n, maxtime, pairwise_interactions, ri, starting_abundances, sigma)

if type(abundances) == int:
    if abundances == -1:
        print("Nonviable")
    else:
        print("Encountered error")
else:
    print("The final abundances are: ")
    print(abundances[-1])
    abundances_line_chart(n, maxtime, abundances)
    #interactions_heatmap(n, pairwise_interactions)

