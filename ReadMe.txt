Simulated data and models to implement the approach presented in the paper under linear scenario for: 
- simulation design 2, with sparse spatial coverege under MAR mechanism; 
- simulation design 3, with sparse spatial coverege under MNAR mechanism.

In particular, the models implement the design stage (Sections 3.2 and 3.3) and the Analysis stage (Sections 3.4).
We provide two examples of code to run the models in R (stored in the folder "DATA_MAR"). 
The code is for the data generated under MAR mechanism, however its use on data generated under 
MNAR mechanism is straightforward to implement.

The simulated data and related shapefiles are available at this link: https://figshare.com/projects/Simulated_data_and_codes_for_the_generalized_EPS/64382

The variables of the simulated data sets are as follows:
	
Variable		Definition
B			Inverse of the correlation matrix used in generating the individual-level variables
C			One ecological level confounder (here available for all the areas)
C.Nomis			One ecological level confounder for the in-sample areas
eps.mean		Posterior mean for the generalized EPS (estimated in the Design stage)
eps.var			Posterior variance for the generalized EPS (estimated in the Design stage)
EPSTrue			True generalized EPS 
EPSTrue.Nomis		True generalized EPS for the in-sample areas
X			Ecological level exposure (here available for all the areas)
X.Nomis			Ecological level exposure for the in-sample areas
Y			Observed number of cases for each area
adjL			Adjacent areas for each area (for entire London - to be used in the ICAR prior in the Analysis stage)
adjSub3			Adjacent areas for each area (for in-sample areas - to be used in the multivariate ICAR prior in the Design stage)
datM			Five mixed-type ecological confounders, which we assume to be unmeasured
datm			Five mixed-type individual-level confounders from survey (here, available for all the areas)
datm.Nomis		Five mixed-type individual-level confounders from survey, available only for the in-sample areas
datm.Nomis.forBugs 	Five-dimension mixed-type individual-level confounders from survey, available only for the in-sample areas for BUGS
expected		Expected number of cases for each area
IndMis			Missing value indicator for individual-level confounders (0=value available; 1=value missing)
nsim			Number of simulated data sets  (here 100)
numareas		Number of areas (for entire London)
numareasNOmis		Number of areas with data on individual-level confounders (i.e. in-sample areas)
numind			Number of subjects for each area with information on individual-level confounders (here 20)
numL			Number of neighbours for each area (for entire London - to be used in the ICAR prior in the Analysis stage)
numSub3			Number of neighbours for each area (for in-sample areas - to be used in the multivariate ICAR prior in the Design stage)
weights3		Weights associated with each pair of areas (for in-sample areas - to be used in the multivariate ICAR prior in the Design stage)
weightsL		Weights associated with each pair of areas (for entire London - to be used in the ICAR prior in the Analysis stage)


