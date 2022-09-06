'''Script for testing out community parameters to be used in simulations.
Simulates one bacterial community with eneralized Lotka-Volterra dynamics. Only pairwise
interactions. Calls appropriate functions for generating interaction values, starting abundances
and for generating trajectories. To modify simulations, modify input values for functions. Standard
deviation input for generating interactions should be in the desired range, but for generating
intrinsic growth rates and pairwise interactions the standard deviation input should be as a fraction
of the given mean. For example if mean is 100 CFU and desired std is 10 CFU the input should be 0.1

Once desired community structure is found, copy community seeds and parameters to data generation.
For easy application the area to copy has been marked with comments "copypaste start" and "copypaste end"
'''

import numpy as np
import pandas as pd
from parameters import generate_growth_rates, generate_interactions, generate_abundances, add_sparcity, adjust_selfinteractions_for_carrying_capacity
from sim_gLV_RK import gLV_RK
from Cyclic_games import random_win_lose_system
from graphics import abundances_line_chart

#starting abundances are defined for community finder only, the repeats will have their own definition
#sa_seed defines how starting abundances are drawn
sa_seed = 55
#starting abundances
sa_mean = 100
sa_std = 0.1

#copypaste start

#Number of species
n = 20
#Time settings for simulation
maxtime = 30
time_increment = 0.1

#Seed values for the random draw of community. Same seeds -> same community
#species_seed chooses how bacteria growth rate and selfinteraction are randomly drawn
species_seed = 23
#interact_seed chooses how interactions are randomly drawn
interact_seed = 7372

#Following parameters determine distributions parameters are drawn from

#Intrinsic growthrate
ri_mean = 0.4
ri_std = 0.1

#Carrying capacities in monoculture, which define selfinteractions
cc_mean = 2000
cc_std = 0.3

#overall interaction strength
interact_mean = 0
interact_std = 0.0004


#Competition settings:

#Add sparcity to interactions, True/False, if True some interactions in the originally drawn pairwise interactions are set to zero
sparcity = False
#sparcity_amount is equal to fraction of interactions set to zero after the initial draw
sparcity_amount = 0.2

#Settings for competition between species
#comp_std defines standard deviation between interactions of winner-loser pairs in competition system, mean is set at winner interaction, loser randomly drawn
comp_std = 0.05
#comp_sparcity defines if competition model is applied to the entire community, sparcity < 1 decribes fraction left out and fraction >= 1 is number of species left out
comp_sparcity=0
#off_target_interaction describes if species interactions left outside the system (comp_sparcity) are set to zero (when False) or kept as is (when True)
off_target_interactions=False

#copypaste end


#Don't touch the following seeds
np.random.seed(species_seed)
ri_seed = np.random.randint(0, 2**31)
cc_seed = np.random.randint(0, 2**31)
np.random.seed(interact_seed)
interactions_seed = np.random.randint(0, 2**31)
interactions_sparcity_seed = np.random.randint(0, 2**31)
competition_interactions_seed = np.random.randint(0, 2**31)
competition_sparcity_seed = np.random.randint(0, 2**31)
competition_magnitude_seed = np.random.randint(0, 2**31)


ri = generate_growth_rates(n, ri_mean, ri_seed, ri_std)

starting_abundances = generate_abundances(n, sa_seed, sa_mean, sa_std)

pairwise_interactions = generate_interactions(n, 1, interactions_seed, interact_mean, interact_std)

if sparcity:
    pairwise_interactions = add_sparcity(pairwise_interactions, sparcity_amount, interactions_sparcity_seed)

pairwise_interactions = random_win_lose_system(n, pairwise_interactions, competition_interactions_seed, comp_sparcity, competition_sparcity_seed, comp_std, competition_magnitude_seed, off_target_interactions)

carrying_capacities = generate_abundances(n, cc_seed, cc_mean, cc_std)

pairwise_interactions = adjust_selfinteractions_for_carrying_capacity(n, pairwise_interactions, ri, carrying_capacities)

abundances = gLV_RK(n, maxtime, time_increment, pairwise_interactions, ri, starting_abundances)

if type(abundances) == int:
    print("Community simulation failed, change parameters")
else:
    if sum(abundances[-1]<10)>0:
        print("One or more communities have very low abundance at endpoint, consider changing parameters.")
    print(abundances[-1])
    abundances_line_chart(n, time_increment, abundances)