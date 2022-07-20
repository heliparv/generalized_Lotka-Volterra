# generalized Lotka-Volterra modelling
Modelling (microbial) communities with generalized Lotka-Volterra dynamics. Two alternate forms of the equation are used, termed here as *simple gLV* and *gLV with carrying capacity* or *gLV with K*

The current simulations are temperamental and especially adding the higher-order interactions results in overflow in most cases.

#### simple gLV
![simple gLV equation](equations/simple_gLV.png)

Where:
- i = species i in group of species
- Ni = species i abundance
- Left side of equation = change of species i abundance per unit of time
- ri = intrinsic growth rate of species i
- Ai<-j = effect of species j on species i
- Bi<-jk = effect of species j and k on species i

#### gLV with K
![gLV with K](equations/gLV_with_K.png)

Where:
- Ki is maximum carrying capacity for species i and it can be found with equation:

![define K](equations/def_K.png) 

- Ktot is maximum total carrying capacity, maximum amount of bacteria the environment supports
- Ntot is the total species abundance

The model scales the species' intrinsic growth rates in relation to the species' carrying capacities so the species growth caps around their individual maximum carrying capacity unless another species affects it positively. The total per capita growth rate for each species is scaled in relation to the total maximum abundance and total maximum carrying capacity so that all species growth slows as the population approaches total maximum carrying capacity.

### Instructions
*index_simple_gLV.py* contains the code needed for initiating a one-off simulation of simple gLV. Input desired parameters to functions and run program.

*index_simple_gLV_HOI.py* contains the code needed for initiating a one-off simulation of simple gLV with higher order interactions. Input desired parameters to functions and run program.

*index_gLVwK.py* contains the code needed for initiating a one-off simulation or gLV with K. Input desired parameters to functions and run program.

*For all simulations*: In random draws the standard deviation values for pairwise interactions should be given in absolute units. For other random draws the standard deviations should be given as fractions of the mean of the distribution used for draw, for example if mean is 100 CFU and wanted std is 10 CFU the input should be 0.1
