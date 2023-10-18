## Data/README
## Purpose: make notes for files in this Data folder

- data.txt: reconstructed individual-level patient data from Hallek et al., by Williams et al.
- 01a_williams2017.txt: organised data from the base case analysis data by Williams et al.
	- basecase15gom1gam2gom3.RData: base case analysis extrapolation results provided by Williams et al.
- 01a_fischer2016.txt: Data generated from /Program/01a_fischer2016.R
	- Fischer2016_Fig1B_RFC.csv: Digitised data for the RFC arm from Fischer et al., Figure 1B
	- Fischer2016_Fig1B_FC.csv: Digitised data for the FC arm from Fischer et al., Figure 1B
	- Fischer2016_Fig2B.json: Digitised data points from Fischer et al., Figure 1B, using the tool from https://apps.automeris.io/wpd/
	- 01a_fischer2016_surv_RFC.txt: Prepare digitised data for the RFC arm from Fischer et al. 
	- 01a_fischer2016_surv_FC.txt: Prepare digitised data for the FC arm from Fischer et al.
	- 01a_fischer2016_nrisk_RFC.txt: Number at risk data for the RFC arm from Fischer et al. 
	- 01a_fischer2016_nrisk_FC.txt: Number at risk data for the FC arm from Fischer et al. 
	- 01a_fischer2016_km_RFC.txt: Kaplan Meier data for the the RFC arm from Fischer et al.
	- 01a_fischer2016_km_FC.txt: Kaplan Meier data for the the FC arm from Fischer et al.
	- 01a_fischer2016_ipd_RFC.txt: Reconstructed individual-level patient for the the RFC arm from Fischer et al.
	- 01a_fischer2016_ipd_FC.txt: Reconstructed individual-level patient for the the FC arm from Fischer et al.
- /popmort/Australia.text, ..., Spain.txt: Population mortality data from 11 countries from the Human Mortality Database
- 01b_popmort.rds: averaged mortality rate weighed by sex across 11 countries
- 04a_rel_survextrap_results_RFC.rds: extrapolation results for the RFC arm by the microsimulation model 
- 04a_rel_survextrap_results_FC.rds: extrapolation results for the FC arm by the microsimulation model
- 07_rel_PSA.RData: Probabilistic sensitivity analysis results from the microsimulation model.
- PSAdata15gom1gam2gom3.RData: Probabilistic sensitivity analysis results provided by Williams et al.

References:
1. Hallek M, Fischer K, Fingerle-Rowson G, Fink A, Busch R, Mayer J, et al. Addition of rituximab to fludarabine and cyclophosphamide in patients with chronic lymphocytic leukaemia: a randomised, open-label, phase 3 trial. The Lancet. 2010 Oct;376(9747):1164–74. 
2. Williams C, Lewsey JD, Briggs AH, Mackay DF. Cost-effectiveness Analysis in R Using a Multi-state Modeling Survival Analysis Framework: A Tutorial. Med Decis Making. 2017 May;37(4):340–52.
3. Fischer K, Bahlo J, Fink AM, Goede V, Herling CD, Cramer P, et al. Long-term remissions after FCR chemoimmunotherapy in previously untreated patients with CLL: updated results of the CLL8 trial. Blood. 2016 Jan 14;127(2):208–15.
4. University of California Berkeley and Max Planck Institute for Demographic Research. The Human Mortality Database [Internet]. Available from: https://www.mortality.org.


