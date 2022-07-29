### LIBRARIES ###
from re import M
import numpy as np
from parameters import generate_growth_rates, generate_interactions, generate_starting_abundances, adjust_selfinteractions, generate_sigma
from graphics import abundances_line_chart
from invasion import multiple_different_species_invasions

### DESCRIPTION ###

"""
Index for performing simulation of multiple mutant invasions based on stochastic version of generalized Lotka-Volterra dynamics. Only pairwise
interactions and the given initial abundances. 
Calls appropriate functions for generating interaction values, starting abundances and for generating trajectories. 
To modify simulations, modify input values for functions. The invasion abundances are assumed to be the initial abundances.
Standard deviation input for generating interactions should be in the desired range, but for generating
intrinsic growth rates, sigma values, and pairwise interactions the standard deviation input should be as a fraction
of the given mean. For example if mean is 100 CFU and desired std is 10 CFU the input should be 0.1
"""

### INITIAL PARAMETERS ###

#Number of species
n = 10
#Maximum simulation time
maxtime = 400


ri = generate_growth_rates(n, 0.4,seed_growth=12, std=0.1)
print("The growth rates: ", ri)
print()

starting_abundances = generate_starting_abundances(n,seed_abundance=12, mean=200, std=0.1)
print("The initial abundances: ", starting_abundances)
print()


pairwise_interactions = generate_interactions(n, 1,seed_interactions=12,seed_sparcity=12,mean= 0,std= 0.0004,sparcity= 0.3)


pairwise_interactions = adjust_selfinteractions(n, pairwise_interactions,seed_selfinter=12 ,mean=-0.0008,std= 0.1)
print("The pairwise interaction grid: \n", pairwise_interactions)
print()

#for total competition dynamics, all negative interactions:
#pairwise_interactions = -1*abs(pairwise_interactions)

sigma = generate_sigma(n,seed_sigma=12, mean=0.1, std=1)
print("The sigma terms of the noise are: \n", sigma)
print()


final_abundances_after_invasions =  multiple_different_species_invasions(n, maxtime, ri, starting_abundances, pairwise_interactions, sigma)
    
print("Final abundances per invasion (rows): \n", final_abundances_after_invasions)
print()
