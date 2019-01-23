This repository contains the R code for the paper  "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity" by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss that has been accepted for publication in Ecology. 

The files and their dependencies are as follows:

Base codes:
===========
base-code.R The base code for simulating the deterministic version of the model. 

stochastic-base-code.R The base code for simulating the stochastic version of the model without overlapping generations. 

stochastic-overlapping-base-code.R The base code for simulating the stochastic version of the model with overlapping generations. 

Figure codes:
=============
FIGURE-1-symmetric.R uses base-code.R to create Figure 1

FIGURE-2-bifs3.R uses base-code.R to create Figure 2

FIGURE-3-whyrho.R uses stochastic-base-code.R to create Figure 3 and Figure 1 in Supplement S2

FIGURE-4-stochastic-bif.R uses stochastic-base-code.R to create Figure 4

FIGURE-5-stochastic-bif-with-storage.R uses stochastic-overlapping-base-code.R to create Figure 5


