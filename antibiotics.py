import numpy as np
from random import sample

"""Functions to add the effect of an antibiotic to the simulation

draw_antibiotic_effect: randomly draws the maximum possible effect of antibiotic
and the Michaelis-Menten constant for each species from normal distributions
with means and standard deviantions as specified in input. Uses given antibiotic
concentration and Michaelis-Menten kinetics to calculate antibiotic effect
at given concentration for the species

add_antibiotic_effect_to_ri: Inputs are number of species, growth rates for species,
inputs needed for draw_antibiotic_effect, and sparcity. Sparcity is either a fraction
or the number of species to drop and the function desides based on that how many
antibiotic effects need to be drawn. Chooses effected species randomly and adds
the effect to the growth rates.

"""

def draw_antibiotic_effect(n, meanmax, stdmax, meanK, stdK, abconc):
    max = -abs(np.random.normal(loc=meanmax, scale=abs(meanmax*stdmax), size=n))
    Km = abs(np.random.normal(loc=meanK, scale=meanK*stdK, size=n))
    eff = np.zeros(n)
    for i in range(0, n):
        eff[i] = max[i]*abconc/(Km[i]+abconc)
    return eff

def add_antibiotic_effect_to_ri(n, ri, meanmax, stdmax, meanK, stdK, abconc, sparcity):
    if sparcity < 1:
        number_effected = int(n*sparcity)
    else:
        number_effected = n-int(sparcity)
    eff = draw_antibiotic_effect(number_effected, meanmax, stdmax, meanK, stdK, abconc)
    effected_species = sample(range(n), number_effected)
    for i in range(number_effected):
        ri[effected_species[i]] += eff[i]
    return ri
