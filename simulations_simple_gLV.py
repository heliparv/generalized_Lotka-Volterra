import numpy as np

""" Functions used to create species trajectories as a function of time with Lotka-Volterra dynamics.

functions:
only_viable_simple_gLV: Generates bacteria abundances as a function of time based on generalized
Lotka-Volterra dynamics. If abundance of some species drops to zero or below the community is
deemed non-viable, the simulation aborts and returns "-1"

simple_gLV_with_extinction: Same as previous but when bacterial abundances run negative, function lets
species go extinct and continues simulation to given maximum time.

test_simple_gLV: No control for negative or zero values, useful for testing how different starting
values affect the calculated abundances.

"""

def only_viable_simple_gLV(n, maxtime, interactions, ri, starting_abundances):
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    unviable = False
    while time < maxtime+1:
        steady = True
        for species in range(0, n):
            change_per_capita = ri[species] + sum(abundances[time-1]*interactions[species])
            if change_per_capita != 0 and steady:
                steady = False
            new_abundance = abundances[time-1][species] + abundances[time-1][species]*change_per_capita
            if new_abundance <= 0:
                unviable = True
                break
            else:
                abundances[time][species] = new_abundance
        if unviable:
            return(-1)
        elif steady:
            break
        else:
            time += 1
    if steady and time != maxtime:
        steady_abundances = np.repeat(np.array([abundances[time]]), maxtime-time, axis=0)
        abundances[time+1:] = steady_abundances
    return abundances

def simple_gLV_with_extinction(n, maxtime, interactions, ri, starting_abundances):
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    while time < maxtime+1:
        i = 0
        steady = True
        while i < len(livespecies):
            change_per_capita = ri[livespecies[i]] + sum(abundances[time-1]*interactions[livespecies[i]])
            if change_per_capita != 0 and steady:
                steady = False
            new_abundance = abundances[time-1][livespecies[i]] + abundances[time-1][livespecies[i]]*change_per_capita
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
        steady_abundances = np.repeat(np.array([abundances[time]]), maxtime-time, axis=0)
        abundances[time+1:] = steady_abundances
    return abundances

def test_simple_gLV(n, maxtime, interactions, ri, starting_abundances):
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    while time < maxtime+1:
        steady = True
        for species in range(0, n):
            change_per_capita = ri[species] + sum(abundances[time-1]*interactions[species])
            if change_per_capita != 0 and steady:
                steady = False
            new_abundance = abundances[time-1][species] + int(np.around((abundances[time-1][species]*change_per_capita), decimals = 9))
            abundances[time][species] = new_abundance
        if steady:
            break
        else:
            time += 1
    if steady and time != maxtime:
        steady_abundances = np.repeat(np.array([abundances[time]]), maxtime-time, axis=0)
        abundances[time+1:] = steady_abundances
    return abundances