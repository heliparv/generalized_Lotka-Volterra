import numpy as np

""" Functions used to create species trajectories as a function of time with Lotka-Volterra dynamics.

functions:
only_viable_trajectories: Generates bacteria abundances as a function of time based on generalized
Lotka-Volterra dynamics. If abundance of some species drops to zero or below the community is
deemed non-viable, the simulation aborts and returns "-1"

#TODO: trajectories: Same as previous but when bacterial abundances run negative, function lets
species go extinct and continues simulation to given maximum time.

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
            change = abundances[species]*change_per_capita
            abundances[time][species] = sum(abundances[species]+change)
            if abundances[time][species] <= 0:
                unviable = True
                break
        if unviable:
            return(-1)
        time += 1
    return abundances

#TODO