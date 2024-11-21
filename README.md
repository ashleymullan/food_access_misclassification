# Linking Potentially Misclassified Healthy Food Access to Diabetes Prevalence
This repository contains `R` code and simulation data to reproduce results from the manuscript by [Mullan, Nguyen, and Lotspeich (2024+)](arxiv-link-will-go-here). 

These simulations rely on the `possum` package, which implements the maximum likelihood approach for covariate misclassification in Poisson regression from the paper. The package can be found in its own repo [here](https://github.com/sarahlotspeich/possum) and installed in `R` as follows:

``` r
# Run once
# install.packages("devtools")
devtools::install_github("sarahlotspeich/possum", ref = "main")
```

## Primary Data Source

The authors gratefully acknowledge [Lotspeich, Mullan, D'Agostino McGowan, and Hepler (2024+)](https://arxiv.org/abs/2405.16385) for the data example used in this paper. For more information on their various data sources, see the documentation in the accompanying [repository](https://github.com/sarahlotspeich/food_access_imputation/blob/main/README.md). The proximities they compute were used to derive indicator variables for binary food access at prespecified radii. For more information, see [our data generation script](https://github.com/ashleymullan/food_access_misclassification/blob/main/Analysis/piedmont_data.R).

## Figures and Tables

### Simulations
[comment]: <> (when I redo these sims, I'll do both sampling strategies, SRS and clever, in the same script)

**Table.** Simulation results under varied two-sided error rates and sampling strategies
- [Script (Run Simulations Locally)](script-location-here)
- [Simulation Results](results-file-location-here)
- [Script (Create Table)](script-location-here)
- [Supplemental Figure (Visualize Results)](script-location-here)

**Table.** Simulation results under varied one-sided error rates and sampling strategies
- [Script (Run Simulations Locally)](script-location-here)
- [Simulation Results](results-file-location-here)
- [Script (Create Table)](script-location-here)
- [Supplemental Figure (Visualize Results)](script-location-here)

**Table.** Simulation results under varied prevalences and sampling strategies
- [Script (Run Simulations Locally)](script-location-here)
- [Simulation Results](results-file-location-here)
- [Script (Create Table)](script-location-here)
- [Supplemental Figure (Visualize Results)](script-location-here)

**Table.** Simulation results under varied prevalence ratios and sampling strategies
- [Script (Run Simulations Locally)](script-location-here)
- [Simulation Results](results-file-location-here)
- [Script (Create Table)](script-location-here)
- [Supplemental Figure (Visualize Results)](script-location-here)

**Table.** Simulation results under varied query proportions and sampling strategies
- [Script (Run Simulations Locally)](script-location-here)
- [Simulation Results](results-file-location-here)
- [Script (Create Table)](script-location-here)
- [Supplemental Figure (Visualize Results)](script-location-here)

### Data Example
**Figure.** Visualize food access at various radii in the Piedmont Triad
- [Script (Generate Figure)](script-location-here)

**Table.** Error rates in the Piedmont Triad data
- [Script (Generate Table)](script-location-here)

**Figure.** Compare estimated diabetes prevalence ratios in the Piedmont Triad under various methods
- [Script (Generate Figure)](script-location-here)

**Figure.** Compare estimated baseline diabetes prevalences in the Piedmont Triad under various methods
- [Script (Generate Figure)](script-location-here)

