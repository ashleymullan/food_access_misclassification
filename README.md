# Linking Potentially Misclassified Healthy Food Access to Diabetes Prevalence
This repository contains `R` code and simulation data to reproduce results from the manuscript by [Mullan, Nguyen, and Lotspeich (2025+)](https://arxiv.org/abs/2505.01465). 

These simulations rely on the `possum` package, which implements the maximum likelihood approach for covariate misclassification in Poisson regression from the paper. The package can be found in its own repo [here](https://github.com/sarahlotspeich/possum) and installed in `R` as follows:

``` r
# Run once
# install.packages("devtools")
devtools::install_github("sarahlotspeich/possum", ref = "main")
```

## Primary Data Source

The authors gratefully acknowledge [Lotspeich, Mullan, D'Agostino McGowan, and Hepler (2025)](https://doi.org/10.1002/sim.70054) for the data example used in this paper. For more information on their various data sources, see the documentation in the accompanying [repository](https://github.com/sarahlotspeich/food_access_imputation/blob/main/README.md). The proximities they compute were used to derive indicator variables for binary food access at prespecified radii. For more information, see [our data generation script](https://github.com/ashleymullan/food_access_misclassification/blob/main/Analysis/case_study.R).

## Figures and Tables

### Simulations

**Table 1.** Simulation results under varied one-sided error rates and sample sizes
- [Script (Run Simulations Locally)](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/one-sided-vary-ppv.R)
- [Simulation Results](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/one-sided-vary-ppv.csv)

**Table 2.** Simulation results under varied query proportions and sample sizes
- [Script (Run Simulations Locally)](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/one-sided-vary-q.R)
- [Simulation Results](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/one-sided-vary-q.csv)

**Table S1.** Simulation results under varied two-sided error rates and sample sizes
- [Script (Run Simulations Locally)](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/two-sided-vary-ppv.R)
- [Simulation Results](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/two-sided-vary-ppv.csv)

**Table S2.** Simulation results under varied prevalences and sample sizes
- [Script (Run Simulations Locally)](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/one-sided-vary-prev.R)
- [Simulation Results](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/one-sided-vary-prev.csv)

**Table S3.** Simulation results under varied prevalence ratios and sample sizes
- [Script (Run Simulations Locally)](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/one-sided-vary-prevrat.R)
- [Simulation Results](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/one-sided-vary-prevrat.csv)

The results from the simulations can be processed into plots with this [script](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/process-sims-into-plots.R) and tables with this [script](https://github.com/ashleymullan/food_access_misclassification/blob/main/Simulations/process-sims-into-tables.R).


### Data Example

This [script](https://github.com/ashleymullan/food_access_misclassification/blob/main/Analysis/case_study.R) generates the version of the data used for analysis in this manuscript and creates Figures 1 (maps of food access in the Piedmont Triad under one-mile and half-mile radii) and 2 (comparison of model results for Piedmont Triad data under four competing analysis methods).


