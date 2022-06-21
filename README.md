# generalized Lotka-Volterra modelling
Modelling (microbial) communities with generalized Lotka-Volterra dynamics.

The current model tends towards trajectories where species abundance dip quickly below zero, so it has to be tweaked.

#### Instructions
*index.py* contains the code needed for initiating a one-off simulation. Input desired parameters to functions and run program.

*parameters.py* contains functions needed to create the initial values for the simulation, such as interactions between species and individual growth rates for species.

*trajectories.py* contains functions needed for carrying out the simulation with the generated initial values.

*graphics.py* contains functions to visualize the generated data.

#### Project roadmap
- Simulation of n species with pairwise interactions
- Third order interactions added to project
- Generalized high-order interactions, order specified in the function call
- Automate creating a dataset of set number of repeats with specified parameters 
