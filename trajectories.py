import numpy as np

""" Functions used to create species trajectories as a function of time with Lotka-Volterra dynamics.

functions:
only_viable_trajectories: Generates bacteria abundances as a function of time based on generalized
Lotka-Volterra dynamics. If abundance of some species drops to zero or below the community is
deemed non-viable, the simulation aborts and returns "-1"

trajectories_with_extinction: Same as previous but when bacterial abundances run negative, function lets
species go extinct and continues simulation to given maximum time.

test_trajectories: No control for negative or zero values, useful for testing how different starting
values affect the calculated abundances.

"""

#TODO separate birth and death rate in model
#TODO take lag times into account

def only_viable_trajectories(n, maxtime, interactions, ri, starting_abundances):
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    unviable = False
    while time < maxtime+1:
        for species in range(0, n):
            change_per_capita = ri[species] + sum(abundances[time-1]*interactions[species])
            new_abundance = abundances[time-1][species] + abundances[time-1][species]*change_per_capita
            if new_abundance <= 0:
                unviable = True
                break
            else:
                abundances[time][species] = new_abundance
        if unviable:
            return(-1)
        time += 1
    return abundances

def trajectories_with_extinction(n, maxtime, interactions, ri, starting_abundances):
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    livespecies = list(range(0,n))
    while time < maxtime+1:
        i = 0
        while i < len(livespecies):
            print(livespecies)
            print(i)
            print(livespecies[i])
            change_per_capita = ri[livespecies[i]] + sum(abundances[time-1]*interactions[livespecies[i]])
            new_abundance = abundances[time-1][livespecies[i]] + abundances[time-1][livespecies[i]]*change_per_capita
            if new_abundance <= 0:
                del livespecies[i]
            else:
                abundances[time][livespecies[i]] = new_abundance
                i+=1
        time += 1
    return abundances

def test_trajectories(n, maxtime, interactions, ri, starting_abundances):
    abundances = np.zeros((maxtime+1, n))
    abundances[0] = starting_abundances
    time = 1
    while time < maxtime+1:
        for species in range(0, n):
            change_per_capita = ri[species] + sum(abundances[time-1]*interactions[species])
            new_abundance = abundances[time-1][species] + abundances[time-1][species]*change_per_capita
            abundances[time][species] = new_abundance
        time += 1
    return abundances