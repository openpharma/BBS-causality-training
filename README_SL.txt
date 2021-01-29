####################################################################################################
####################################################################################################
## Accompying R code for the BBS training "A gentle introduction into causal thinking", Feb 2 2021
## Authors: Sabine Lauer, Shengchun Kong, Dominik Heinzmann* (* dominik.heinzmann@roche.com)  
## Date: Jan 29, 2021
####################################################################################################
####################################################################################################


FILES
1) example_data.csv: Example data
2) ADA_seminar_example_v2.R: Is the main file which needs to be executed
3) ADA_functions_v2.R: Contains all relevant help functions


OUTPUT FILES
A) report_Example study_OS.csv: Contains the results of the principal stratum approach, as well as unadjusted results (ie ADA groups vs entire control) and ADA+ vs ADA- (only treatment arm)
B) KM_ADAneg.png & KM_ADApos.png: Kaplan-Meier curves for the ADA negative respectively positive group versus the appropriate (re-weighted) control
C) KM_Controlgrps.png: Kaplan-Meier curves of entire control and the re-weighted (principal stratum) controls for both ADA groups
D) KM_ADAgrps.png: Kaplan-Meier curves of ADA groups
E) samplesize_Example study_OS.csv: Contains sample sizes for the ADA groups (landmark groups), assuming that there are no missing covariates
F) Balance_negative.png & Balance_positive.png: Balance diagnostics plots for both ADA groups
G) Feature_count.csv: Results of supplemental model diagnostics as introduced in the course 
H) warningnum_Example study_OS.csv: Contains the results of the bootstrap run. If numbers of warnings is high, one needs to check the bootstrap samples - eg some covariate categories with low incidence might often not be sampled from one of the ADA groups. If needed, list of covariates needs some considerations and adjustments. 



SOME GUIDANCE HOW TO ADJUST FOR YOUR OWN EXAMPLE
a) Structure your example similar to "example_data.csv"; note that the ADA status "ADALM" already contains the landmark ADA status and that the population has already been reduced to patients still alive and not censored at the landmark (OSLM="Y"). So you will need to apply the landmark approach before (or add corresponding code at the beginning of ADA_seminar_example_v2.R).  
b) Check "ADA_seminar_example_v2.R" to adjust covariate names in the first part of the code (wherever used, eg list of covariates, model formulas, and dichotomization for balance diagnostics) 
c) Check "ADA_seminar_example_v2.R" to adjust for the naming of the categories in ADALM (in case of a landmark different from week 4)
d) Change the directory to one suitable for you, see function setwd() in "ADA_seminar_example_v2.R"
d) Execute "ADA_seminar_example_v2.R" after changes per a)-d) 
    
