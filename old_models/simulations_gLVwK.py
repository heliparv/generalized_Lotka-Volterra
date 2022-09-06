import numpy as np

""" Functions used to create species trajectories as a function of time with Lotka-Volterra dynamics.
All functions test for achieving steady state and if there's no change to any of the abundances during
a simulation round the simulation is interrupted and steady state abundances are filled in for the
remainder of the simulation time.

functions:
only_viable_gLVwK: Generates bacteria abundances as a function of time based on generalized
Lotka-Volterra dynamics. If abundance of some species drops to zero or below the community is
deemed non-viable, the simulation aborts and returns "-1". If overflow encountered returns "-2"
in stead of abundances.

gLVwK_with_extinction: Same as previous but when bacterial abundances run negative, function lets
species go extinct and continues simulation to given maximum time. If overflow encountered returns "-2"
in stead of abundances.

test_gLVwK: No control for negative or zero values, useful for testing how different starting
values affect the calculated abundances. If overflow encountered, prints out values causing overflow
and returns abundances calculated so far.

calculate_abundance_product_vector: calculates the product of abundances according to the call vector
created for modelling higher order interactions in the model

"""

def only_viable_gLVwK(n, maxtime, time_increment, interactions, ri, carrying_capacities, starting_abundances, total_carrying_capacity):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((max_increments+1, n))
    abundances[0] = starting_abundances
    time = 1
    nonviable = False
    while time < max_increments+1:
        steady = True
        total_carrying_limit = (total_carrying_capacity-sum(abundances[time-1]))/total_carrying_capacity
        for species in range(0, n):
            try:
                change_per_capita = (ri[species]*((carrying_capacities[species]-abundances[time-1][species])/carrying_capacities[species])+ (sum(abundances[time-1]*interactions[species])))*time_increment
                change_per_capita = total_carrying_limit*change_per_capita
                if abs(change_per_capita) > 0.0001 and steady:
                    steady = False
                new_abundance = np.around((abundances[time-1][species] + abundances[time-1][species]*change_per_capita), decimals=6)
            except:
                return -2
            if new_abundance <= 0:
                nonviable = True
                break
            else:
                abundances[time][species] = new_abundance
        if nonviable:
            return -1
        elif steady:
            break
        else:
            time += 1
    if steady and time != max_increments:
        steady_abundances = np.repeat(np.array([abundances[time]]), max_increments-time, axis=0)
        abundances[time+1:] = steady_abundances
    return abundances

def gLVwK_with_extinction(n, maxtime, time_increment, interactions, ri, carrying_capacities, starting_abundances, total_carrying_capacity):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((max_increments+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    while time < max_increments+1:
        i = 0
        steady = True
        total_carrying_limit = (total_carrying_capacity-sum(abundances[time-1]))/total_carrying_capacity
        while i < len(livespecies):
            try:
                change_per_capita = (ri[livespecies[i]]*((carrying_capacities[livespecies[i]]-abundances[time-1][livespecies[i]])/carrying_capacities[livespecies[i]])+ (sum(abundances[time-1]*interactions[livespecies[i]])))*time_increment
                change_per_capita = total_carrying_limit*change_per_capita
                if abs(change_per_capita) > 0.0001 and steady:
                    steady = False
                new_abundance = np.around((abundances[time-1][livespecies[i]] + abundances[time-1][livespecies[i]]*change_per_capita), decimals=6)
            except:
                return -2
            if new_abundance <= 0:
                del livespecies[i]
            else:
                abundances[time][livespecies[i]] = new_abundance
                i+=1
        if steady:
            break
        else:
            time += 1
    if steady and time < max_increments:
        steady_abundances = np.repeat(np.array([abundances[time]]), max_increments-time, axis=0)
        abundances[time+1:] = steady_abundances
    return abundances

def test_gLVwK(n, maxtime, time_increment, interactions, ri, carrying_capacities, starting_abundances, total_carrying_capacity):
    np.seterr(all='raise')
    max_increments = int(maxtime*(1/time_increment))
    abundances = np.zeros((max_increments+1, n))
    abundances[0] = starting_abundances
    time = 1
    while time < max_increments+1:
        steady = True
        total_carrying_limit = (total_carrying_capacity-sum(abundances[time-1]))/total_carrying_capacity
        for species in range(0, n):
            try:
                change_per_capita = (ri[species] *((carrying_capacities[species]-abundances[time-1][species])/carrying_capacities[species])+ (sum(abundances[time-1]*interactions[species])))*time_increment
                change_per_capita = total_carrying_limit*change_per_capita
            except:
                print(f"Time: {time}")
                print(f"species: {species}")
                print(f"abundances t-1: {abundances[time-1]}")
                print(f"interactions: {interactions[species]}")
                print(f"Carrying capasity: {carrying_capacities[species]}")
                print(f"total carrying limit: {total_carrying_limit}")
                return abundances
            if abs(change_per_capita) > 0.0001 and steady:
                steady = False
            try:
                new_abundance = np.around((abundances[time-1][species] + abundances[time-1][species]*change_per_capita), decimals=6)
                abundances[time][species] = new_abundance
            except:
                print(f"Time: {time}")
                print(f"Species: {species}")
                print(f"Abundance t-1: {abundances[time-1][species]}")
                print(f"Change per capita: {change_per_capita}")
                return abundances
        if steady:
            break
        else:
            time += 1
    if steady and time < max_increments:
        steady_abundances = np.repeat(np.array([abundances[time]]), max_increments-time, axis=0)
        abundances[time+1:] = steady_abundances
    return abundances

def calculate_abundance_product_vector(call, abundances):
    product = np.ones(len(call))
    for i in range(0, len(call)):
        for j in call[i]:
            product[i] = product[i]*abundances[j]
    return product