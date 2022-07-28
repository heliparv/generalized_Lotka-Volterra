import numpy as np

""" Functions used to create species trajectories as a function of time with a generalized
Lotka-Volterra model with nutrient dynamics. Change rate of a bacteria species is determined
by it's intrinsic maximum rate of growth, amount of nutrient present, and the effect of self-
interactions and interactions with other species. Only pairwise interactions. Functions test
for achieving steady state and if there's no change to any of the abundances during a round
the simulation is interrupted and steady state abundances are filled in for the remainder of
the simulation time. Amount of nutrient is modelled as a chemostat with fixed inflow of nutrient
and with bacteria abundance dependant loss of nutrient.

functions:
only_viable_nd: Generates bacteria abundances as a function of time based on generalized
Lotka-Volterra dynamics. If abundance of some species drops to zero or below the community is
deemed non-viable, the simulation halts and returns "-1". If overflow encountered returns "-2"
in stead of abundances.

nd_with_extinction: Same as previous but when bacteria abundances run negative, function lets
species go extinct and continues simulation to given maximum time. If overflow encountered returns
"-2" in stead of abundances.

"""

def only_viable_nd(n, maxtime, interactions, ri, Ki, starting_abundances, starting_nutrient, gamma, influx):
    np.seterr(all='raise')
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    nutrient = np.zeros(maxtime+1)
    nutrient[0] = starting_nutrient
    time = 1
    while time < maxtime+1:
        steady = True
        try:
            growth_rates = (ri*nutrient[time-1])/(Ki+nutrient[time-1])
            nutrient[time] = nutrient[time-1] + influx - sum(gamma*growth_rates*abundances[time-1])
            if nutrient[time] != nutrient[time-1]:
                steady = False
        except:
            return -2, -2
        for species in range(0, n):
            try:
                change_per_capita = growth_rates[species] + sum(abundances[time-1]*interactions[species])
                if change_per_capita != 0 and steady:
                    steady = False
                new_abundance = abundances[time-1][species] + int(abundances[time-1][species]*change_per_capita)
            except:
                return -2, -2
            if new_abundance <= 0:
                return -1, -1
            else:
                abundances[time][species] = new_abundance
        if steady:
            break
        else:
            time += 1
    if steady and time != maxtime:
        return add_steady_state_tail(abundances, nutrient, time, maxtime)
    return abundances, nutrient

def nd_with_extinction(n, maxtime, interactions, ri, Ki, starting_abundances, starting_nutrient, gamma, influx):
    np.seterr(all='raise')
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    nutrient = np.zeros(maxtime+1)
    nutrient[0] = starting_nutrient
    time = 1
    livespecies = list(range(0,n))
    while time < maxtime+1:
        i = 0
        steady = True
        try:
            growth_rates = (ri*nutrient[time-1])/(Ki+nutrient[time-1])
            nutrient[time] = nutrient[time-1] + influx - sum(gamma*growth_rates*abundances[time-1])
            if nutrient[time] != nutrient[time-1]:
                steady = False
        except:
            return -2, -2
        while i < len(livespecies):
            try:
                change_per_capita = growth_rates[livespecies[i]] + sum(abundances[time-1]*interactions[livespecies[i]])
                if change_per_capita != 0 and steady:
                    steady = False
                new_abundance = abundances[time-1][livespecies[i]] + int(abundances[time-1][livespecies[i]]*change_per_capita)
            except:
                return -2, -2
            if new_abundance <= 0:
                del livespecies[i]
            else:
                abundances[time][livespecies[i]] = new_abundance
                i+=1
        if steady:
            break
        else:
            time += 1
    if steady and time != maxtime:
        return add_steady_state_tail(abundances, nutrient, time, maxtime)
    return abundances, nutrient

def add_steady_state_tail(abundances, nutrient, time, maxtime):
    steady_abundances = np.repeat(np.array([abundances[time]]), maxtime-time, axis=0)
    abundances[time+1:] = steady_abundances
    steady_nutrient = np.repeat(np.array([nutrient[time]]), maxtime-time)
    nutrient[time+1:] = steady_nutrient
    return abundances, nutrient
