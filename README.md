Supplementary material C - data package and model code
================
R.M. Link, S. Fuchs, D. Arias Aguilar, C. Leuschner, M. Castillo Ugalde,
J.C. Valverde Otarola, B. Schuldt

## Description

The present R project is part of the digital supplementary material of
Link et al. (2020) – *Tree height predicts the shape of radial sap flow
profiles of Costa-Rican tropical dry forest tree species* – and contains
the Stan code for the three different hierarchical models considered in
the paper (cf. Supplementary Material A for details about the model
structure). The article is available under the following link:
<https://doi.org/10.1016/j.agrformet.2020.107913>.

## Structure

The Stan model code is stored in the `/stan` folder, while the code for
the simple frequentist models used to generate starting values is stored
in `1_HFD_analysis_nlme.R` and the code used to fit the models is found
in `1_HFD_analysis_stan.R`.

The full project structure is as follows:

``` text
/               top-level scripts:
                1_HFD_analysis_nlme.R   -- preliminary nlme-based models for starting values
                2_HFD_analysis_stan.R   -- final model fitting with Stan
/data           dataset used for model fitting
/output         saved model objects, model output etc.
/R              scripts that are sourced by the top level scripts
                /R/get_starting_values_ABC.R   -- functions for randomized starting values based
                                                  on the maximum likelihood estimates
                /R/utility_functions.R         -- utility functions used in the top level scripts                                  
/stan           model code for the Stan models
                /stan/model_A*       -- simple model without species level random effects and 
                                        parameter regressions
                /stan_model_B*       -- model with full random effects but without parm. regressions
                /stan_model_C*       -- full model
```
