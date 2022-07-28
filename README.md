# cosm_parameter_mcmc

uses mcmc to reconstruct posterior of alpha dgp cosmological model fit to supernovae observations

1.: set STEPS in mcmcMethods.h to desired number of steps (default= 1e6)
2.: compile and run with $ make run
3.: resulting chain will be in ./out/chain.txt

chain.txt has three coloumns: 1. omega_m; 2. H0; 3. alpha
