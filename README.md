# Linking Potentially Misclassified Healthy Food Access to Diabetes Prevalence
This repository contains `R` code and simulation data to reproduce results from the manuscript by [Mullan, Nguyen, and Lotspeich (2024+)](arxiv-link-will-go-here). 

These simulations rely on the `possum` package, which implements the maximum likelihood approach for covariate misclassification in Poisson regression from the paper. The package can be found in its own repo [here](https://github.com/sarahlotspeich/possum) and installed in `R` as follows:

``` r
# Run once
# install.packages("devtools")
devtools::install_github("sarahlotspeich/possum", ref = "main")
```

## Primary Data Source (Analysis)

The authors gratefully acknowledge [Lotspeich, Mullan, D'Agostino McGowan, and Hepler (2024+)](https://arxiv.org/abs/2405.16385) for the data example used in this paper. For more information on how it was put together, see the documentation in the accompanying [repository](https://github.com/sarahlotspeich/food_access_imputation/blob/main/README.md). The proximities they compute were used to derive indicator variables for binary access at prespecified radii.
