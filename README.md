# Power-Lindley Generalized Pareto (PL-GPD) distribution

R functions for the PL-GPD composite distribution, from my M.Sc. Statistics thesis,
*Parameter Estimation of a New Family of Distribution for Heavy Tailed Dataset with Application*
(University of Lagos, supervisor: Dr E. B. Nkemnole).

## The model
Data such as financial returns mix many ordinary values with a few extreme ones. A composite
distribution uses one distribution for the bulk and another for the tail, joined at a threshold `u`.
Here the bulk is a Power Lindley distribution and the tail is a Generalized Pareto distribution (GPD).

## What is in this repository
`plgpd_functions.R` contains:
- density, distribution and quantile functions of the Power Lindley and Generalized Pareto distributions
- density, distribution and quantile functions of the PL-GPD (`dPLGPD`, `pPLGPD`, `qPLGPD`)
- example plots of the PL-GPD density and distribution function for several parameter settings

Parameters: `u` threshold, `a` and `b` Power Lindley parameters (alpha, beta), `r` GPD shape (gamma), `c` GPD scale.
The mixing weight `w` is computed inside each function. In these functions the scale `c` is an input.

## Not included
The maximum likelihood fitting and the S&P 500 log-returns analysis reported in the thesis are not in this repository.

## Run it
In R, run: `source("plgpd_functions.R")`
Packages loaded by the script: LambertW, lamW, and others listed at the top of the file.

## Publication
E. B. Nkemnole and N. G. Okoacha (2025). Power-Lindley Generalized Pareto Distribution: A New Approach
for Modeling Heavy-Tailed Data. *Research and Analysis Journals*, 8(04), 06-20.
doi:10.18535/raj.v8i04.514

## Author
Nkechi Grace Okoacha, Pan-Atlantic University, Lagos.
