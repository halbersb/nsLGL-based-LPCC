# LPCC_based_nsLGL
LPCC-based non-stationary local-to-global learning
implemented in Matlab.

This package includes all the relevant dependencies

To run the code:
1. download this folder
   (a python script to clone from "anonymous.4open.science" can be found here https://github.com/ShoufaChen/clone-anonymous4open)
2. open in matlab (2015 or later)
3. run main.m

data should be in the format: (S,O*T)
where S in the number of sequences, O is the number of observed variables, and T is the number of time-slices

**To run the "latent variable time-varying graphical lasso" (LTGL) go to git "fdtomasi/regain"

running example.xlsx demonstrates how LPCC-based nsLGL algorithm transforms local graphs into a global nsDBN