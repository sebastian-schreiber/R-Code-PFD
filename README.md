This repository contains the R code for the paper  "When rarity has costs: coexistence under positive frequency-dependence and environmental stochasticity" by  Sebastian J. Schreiber, Masato Yamamichi and Sharon Y. Strauss that has been accepted for publication in Ecology. 

The files and their dependencies are as follows:

Base codes:
===========
base-code.R The base code for simulating the deterministic version of the model. 

base-code-stochastic.RThe base code for simulating the stochastic version of the model without overlapping generations. 

base-code-stochastic-overlapping.R The base code for simulating the stochastic version of the model with overlapping generations. 

Figure codes:
=============
FIGURE-1-symmetric.R uses base-code.R to create Figure 1

FIGURE-2-bifs3.R uses base-code.R to create Figure 2

FIGURE-3-why-rho.Ruses stochastic-base-code.R to create Figure 3 and Figures S-1A,B in Supplement S3

FIGURE-4-stochastic-bif.R uses stochastic-base-code.R to create Figure 4

FIGURE-5-stochastic-bif-with-storage.R uses stochastic-overlapping-base-code.R to create Figure 5

FIGURE-6-stochasticbifs4.R uses stochastic-base-code.R to create Figure 6

FIGURE-S3-1C.R uses stochastic-overlapping-base-code.R to create Figure S-1C in Supplement S3

