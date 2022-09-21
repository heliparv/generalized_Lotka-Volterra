''' Index for data creation for a community of microbes with a generalized Lotka-Volterra model.
Bacteria community can be built in community_finder.py by altering input values and seeds,
once community has been found it can be copypasted here, for ease of use comments mark
"copypaste start" and "copypaste end" to show which values should be copied to achieve the same
community structure.

Simulates a community with all species present and leave-one-out communities.

Output settings are before copypaste area.

The simulation does not test for failed communities, so even in in some simulations the abundance
of one species or more drops to zero the data is still collected.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from parameters import generate_growth_rates, generate_interactions, generate_abundances, add_sparcity, adjust_selfinteractions_for_carrying_capacity, generate_loo_starting_abundances
from sim_gLV_RK import gLV_RK
from Cyclic_games import random_win_lose_system, random_competition_system
from graphics import abundances_line_chart
from create_data_outputs import initiate_df, output, relative_abundances, addable_to_frame

#Output settings
#choose number of sampled time points per experiment
no_of_samples = 23
#In addition to numerical data the simulations can save figures of the community simulations.
# These include the whole simulation, most of which is not available in the numerical data
save_figures_of_communities = True
#name of the output used for files
output_name = "output"

#leave-one-out simulations draw the starting abundances from the same distribution for all experiments
#sa_seed defines how starting abundances are drawn
sa_seed = 55
#starting abundances mean and standard deviation
sa_mean = 100
sa_std = 0.05

#copypaste start

#Number of species
n = 23
#Time settings for simulation
maxtime = 25
time_increment = 0.05

#Seed values for the random draw of community. Same seeds -> same community
#species_seed chooses how bacteria growth rate and selfinteraction are randomly drawn
species_seed = 567
#interact_seed chooses how interactions are randomly drawn
interact_seed = 73112

#Following parameters determine distributions parameters are drawn from

#Intrinsic growthrate
ri_mean = 0.6
ri_std = 0.2

#Carrying capacities in monoculture, which define selfinteractions
cc_mean = 2000
cc_std = 0.2

#overall interaction strength
interact_mean = -0.00003
interact_std = 0.000006


#Competition settings:

#Add sparcity to interactions, True/False, if True some interactions in the originally drawn pairwise interactions are set to zero
sparcity = False
#sparcity_amount is equal to fraction of interactions set to zero after the initial draw
sparcity_amount = 0

#Settings for competition between species
#comp_std defines standard deviation between interactions of winner-loser pairs in competition system, mean is set at winner interaction, loser randomly drawn
comp_std = 0.05
#comp_sparcity defines if competition model is applied to the entire community, sparcity < 1 decribes fraction left out
comp_sparcity=0.7
#off_target_interaction describes if species interactions left outside the system (comp_sparcity) are set to zero (when False) or kept as is (when True)
off_target_interactions=False

#copypaste end
np.random.seed(sa_seed)
whole_com_sa_seed = np.random.randint(0, 2**31)
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

whole_com_starting_abundances = generate_abundances(n, whole_com_sa_seed, sa_mean, sa_std)

starting_abundances = generate_loo_starting_abundances(n, sa_seed, sa_mean, sa_std)

pairwise_interactions = generate_interactions(n, 1, interactions_seed, interact_mean, interact_std)

if sparcity:
    pairwise_interactions = add_sparcity(pairwise_interactions, sparcity_amount, interactions_sparcity_seed)

pairwise_interactions = random_competition_system(n, pairwise_interactions, comp_sparcity, competition_sparcity_seed, comp_std, competition_magnitude_seed, off_target_interactions)

carrying_capacities = generate_abundances(n, cc_seed, cc_mean, cc_std)

pairwise_interactions = adjust_selfinteractions_for_carrying_capacity(n, pairwise_interactions, ri, carrying_capacities)

param_names_column = ["sa_seed", "sa_mean", "sa_std", "n", "maxtime", "time_increment", "species_seed", "interact_seed", "ri_mean", "ri_std", "cc_mean", "cc_std", "interact_mean", "interact_std", "sparcity", "sparcity_amount", "comp_std", "comp_sparcity", "off_target_interactions"]

param_column = [sa_seed, sa_mean, sa_std, n, maxtime, time_increment, species_seed, interact_seed, ri_mean, ri_std, cc_mean, cc_std, interact_mean, interact_std, sparcity, sparcity_amount, comp_std, comp_sparcity, off_target_interactions]

param_df = pd.DataFrame(columns=["parameter name", "parameter value"])
param_df["parameter name"] = param_names_column
param_df["parameter value"] = param_column

columns, absolute_abund_df = initiate_df(n)
columns, relative_abund_df = initiate_df(n)

sampled_times = np.linspace(0, int(maxtime*(1/time_increment)), num=no_of_samples).astype(int)

abundances = gLV_RK(n, maxtime, time_increment, pairwise_interactions, ri, whole_com_starting_abundances)

if save_figures_of_communities:
    plt.plot(list(range(int(maxtime*(1/time_increment))+1)), abundances)
    plt.title(f"Whole community")
    plt.xlabel("Time")
    plt.ylabel("Abundance")
    plt.savefig(f"{output_name}_whole_community.png")
    plt.clf()

sample = np.take(abundances, sampled_times, axis=0)
absolute_abund_df = absolute_abund_df.append(addable_to_frame(columns, 0, sampled_times, sample))
rel_data = relative_abundances(sample)
relative_abund_df = relative_abund_df.append(addable_to_frame(columns, 0, sampled_times, rel_data))

for i in range(len(starting_abundances)):
    abundances = gLV_RK(n, maxtime, time_increment, pairwise_interactions, ri, starting_abundances[i])
    if save_figures_of_communities:
        plt.plot(list(range(int(maxtime*(1/time_increment))+1)), abundances)
        plt.title(f"Species left out: {i+1}")
        plt.xlabel("Time")
        plt.ylabel("Abundance")
        plt.savefig(f"{output_name}_leave_{i+1}_out.png")
        plt.clf()
    sample = np.take(abundances, sampled_times, axis=0)
    absolute_abund_df = absolute_abund_df.append(addable_to_frame(columns, i+1, sampled_times, sample))
    rel_data = relative_abundances(sample)
    relative_abund_df = relative_abund_df.append(addable_to_frame(columns, i+1, sampled_times, rel_data))

output(output_name, param_df, absolute_abund_df, relative_abund_df, n, pairwise_interactions)
