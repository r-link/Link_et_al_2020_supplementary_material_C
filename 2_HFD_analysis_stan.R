################################################################################
#	
# HFD data analysis with rstan   
#
# Aim: generate and inspect HMC samples from non-linear hierarchical 
#      model for vertical sap flow profiles
#
# Author: Roman Link    Contact: roman.link@plant-ecology.de
#
################################################################################

# Sections:

# 1. Data loading and preparation
# 2. Fit model A
# 3. Fit model B
# 4. Fit model C
# 5. Inspect model output
# 6. Inspect variance decomposition
# 7. Plot expected profiles and confidence intervals

# load packages ----------------------------------------------------------------
# list of packages 
pkgs <-c("tidyverse", # framework for data manipulation
         "magrittr",  # additional piping operators
         "nlme",      # nonlinear mixed effects models (needed for starting value script)
         "rstan",     # stan interface for R (make sure to have stan installed on your system!)
         "shinystan"  # shiny model output for rstan
        )        
# check for existence of packages and install if necessary
to_install<-pkgs[!(pkgs %in% installed.packages()[,1])]
if (length(to_install)>0)  for (i in seq(to_install)) install.packages(to_install[i])
# load all required packages 
for (i in seq(pkgs)) require(pkgs[i], character.only = T)

# set rstan option settings
rstan_options(auto_write = TRUE)              # automatically save compiled c++ model objects in same folder as .stan files
options(mc.cores = parallel::detectCores())   # use as many parallel cores as supported by the system


# 1. Data loading and preparation ----------------------------------------------
# load script with helper functions
source("R/utility_functions.R")

# load sap flow dataset
set.seed(101) # set seed for plotting order
data <- read_csv("data/radial_profile_data_clean.csv") %>%
  mutate(species  = factor(species),
         tree1 = factor(tree, levels = sample(unique(tree))), # randomly reorder species variable (for plotting)
         tree     = factor(tree),
         WD       = scale1(WD), # scale1 is a wrapper around scale that drops attributes (see utility functions)
         height   = scale1(height),
         ASI      = scale1(ASI),
         treedate = factor(paste(species, tree, date, sep = "/"))) 


# define list with all data needed to fit the models in /stan/*
data_fit <- data %$%
  list(I        = nrow(data),           # number of observations
       flux     = flux,                 # daily averages of sap flux rates
       depth    = reldepth,             # sensor installation depths
       J        = nlevels(treedate),    # number of tree:date combinations 
       K        = nlevels(tree),        # number of trees
       L        = nlevels(species),     # number of species
       tree     = as.numeric(tree),     # tree ID
       treedate = as.numeric(treedate), # measurement days for each tree
       spec     = as.numeric(species),  # species ID
       unique   = which(!duplicated(tree)), # first values of each tree (for variance decomposition)
       X        =  model.matrix(~ WD + height + ASI)  # predictor matrix for the fixed effects
  )

# load script with function for randomized starting values
source("R/get_starting_values_ABC.R")


# 2. Fit model A --------------------------------------------------------------
# BE CAREFUL
# Model fitting takes a lot of time, and the exported model objects consume a lot 
# of hard disk space. Run the following lines at own risk!

# fit model with rstan
seed_A <- 20168 # predefined seed for reproducibility
stan_fit_A <- stan(
                 # path to stan model code
                 file = "stan/model_A_spec_only.stan", 
                 # list with data used for modelling
                 data   = data_fit, 
                 # number of iterations per chain
                 iter   = 10000, 
                 # number of iterations discarded for warmup 
                 warmup =  5000,
                 # number of chains
                 chains =     4,
                 # seed used for model fitting
                 seed   = seed_A,
                 # function for starting values
                 init   = startfun_A,
                 # untracked (temporary) parameters
                 pars   = c("u_raw", "v_raw", "z_spec", "z_spec_1", "L_Omega_spec"),
                 # flag that declares that parameters named in "pars" have to be discarded
                 include = FALSE,
                 # print progress?
                 verbose = TRUE,
                 # control arguments for the sampling process 
                 control = list(adapt_delta    = 0.999,
                                max_treedepth  = 15,
                                stepsize       = 0.01),
                 # how often shall the progress be printed?
                 refresh = 10)

# export results 
save(stan_fit_A, file = "output/stan_fit_A.RData")


# 3. Fit model B --------------------------------------------------------------
seed_B <- 924805 # predefined seed for reproducibility
stan_fit_B <- stan(file = "stan/model_B_tree_RE.stan", 
                 data   = data_fit,
                 iter   = 10000, 
                 warmup =  5000,
                 chains =     4,
                 seed   = seed_B,
                 init   = startfun_B,
                 pars   = c("u_raw", "z_tree", "z_spec","z_tree_1", "z_spec_1",
                            "L_Omega_tree", "L_Omega_spec"),
                 include = FALSE,
                 verbose = TRUE,
                 control = list(adapt_delta    = 0.999,
                                max_treedepth  = 15,
                                stepsize       = 0.01),
                 refresh = 10)

# export results
save(stan_fit_B, file = "output/stan_fit_B.RData")


# 4. Fit model C --------------------------------------------------------------
# fit model with rstan
seed_C <- 658739 # predefined seed for reproducibility
stan_fit_C <- stan(file = "stan/model_C_tree_RE_and_fixed_effects.stan", 
                 data   = data_fit,
                 iter   = 10000, 
                 warmup =  5000,
                 chains =     4,
                 seed   = seed_C,
                 init   = startfun_C,
                 pars   = c("u_raw", "z_tree", "z_spec","z_tree_1", "z_spec_1",
                            "L_Omega_tree", "L_Omega_spec"),
                 include = FALSE,
                      append_samples = FALSE, 
                 verbose = TRUE,
                 control = list(adapt_delta    = 0.999,
                                max_treedepth  = 15,
                                stepsize       = 0.01),
                 refresh = 10)


# export results
save(stan_fit_C, file = "output/stan_fit_C.RData")


# 5. Inspect model output ------------------------------------------------------
# the following section is identical for all models - to avoid redundancy,
# just specify the model to inspect here:
stan_fit <- stan_fit_C 

# shiny stan model output (BE CAREFUL - MAY TAKE VERY LONG)
launch_shinystan(stan_fit)

# summary of most important parameters
# for model A and B (uncomment if necessary)
# parnames <- c("mu0", "k0", "mult0", "sigma", "tau_u", "tau_v", "tau_w",
#               "Omega_spec", "hier_scale")

# for model C
parnames <- c("beta_mu", "beta_k", "mult0", "sigma", "tau_u", "tau_v", "tau_w",
              "Omega_tree", "Omega_spec", "hier_scale")

# get summaries for the most important parameters
(mcmcsumm <- summary(stan_fit, pars = parnames)$summary %>% as.data.frame() %>% .[.$mean != 1,])

# trace plot of most important parameters
traceplot(stan_fit, inc_warmup = F, parnames)

# plot estimated parameters
pars <- extract_chain(stan_fit, parnames[!grepl("Omega", parnames, perl = T)])
pars %>%
  gather("parameter", "estimate") %>%
  mutate(class = ifelse(grepl("beta|k0|mu0|mult0", parameter, perl = T), "Parameter", "Variance")) %>%
  ggplot(aes(x = parameter, y = estimate)) + 
  geom_violin(alpha = 0.3, fill = "steelblue", scale = "width")+
  geom_hline(yintercept = 0, lty = 2)+
  stat_summary(geom = "pointrange", fun.data = mci) + 
  facet_wrap(~class, scales = "free")

# compare parameter estimates of fixed effects to ML estimates (only model C)
(meanpars <- select(pars, mult0,contains("beta")) %>% 
    mutate(mult0 = mult0 - log(sd(data_fit$flux))) %>%
  sapply(mean))
fixef(model_C)
plot(meanpars ~ fixef(model_C)); abline(0, 1, col = "grey")
# --> estimated parameters are very similar to the maximum likelihood estimates
#     but the confidence bounds are much wider (almost everything was "significant"
#     in the nlme version of the model due to problems with the number of df for testing)


# 6. Inspect variance decomposition  --------------------------------------------
# variance components
# a) mu
(mu_vars <- extract_chain(stan_fit, "mu_vardecomp") %>%  
   colMeans %>%
   set_names(c("% var. explained by fixed effects", 
               "% var. explained by species RE",
               "% var. explained by tree RE"))) 
# b) disp
(k_vars  <- extract_chain(stan_fit, "k_vardecomp") %>% 
    colMeans%>%
    set_names(c("% var. explained by fixed effects", 
                "% var. explained by species RE",
                "% var. explained by tree RE"))) 
# c) pseudo-rsq
(pseudo  <- extract_chain(stan_fit, "y_vardecomp") %>% 
    colMeans %>%
    set_names(c("% var. explained by tree level preds.", 
                "% var. explained by day RE",
                "% unexplained variance"))) 


# 7. Plot expected profiles and confidence intervals----------------------------
# ...preparation ------------
# extract day-level profiles
EP <- extract_chain(stan_fit, "exp_profile") * sd(data$flux) 
# calculate summary statistics
EP_summary <- t(apply(EP, 2, function(x)
  c(mean = mean(x), quantile(x, c(0.025, 0.975)))))
colnames(EP_summary) <- c("mean","ymin","ymax")

# get average profiles
EP_tree <- EP / exp(extract_chain(stan_fit, "u")[,data_fit$treedate])
# calculate summary statistics
EP_tree_summary <- t(apply(EP_tree, 2, function(x)
  c(mean = mean(x), quantile(x, c(0.025, 0.975)))))
colnames(EP_tree_summary) <- c("mean_tree","ymin_tree","ymax_tree")

# bind to raw data
datap <- cbind(data, EP_summary, EP_tree_summary) %>% as_tibble %>%
  mutate(tree = forcats::fct_reorder(tree, as.numeric((species))))

# ...species level -------
ggplot(datap)+
  geom_point(aes(x = depth, y = flux, col = tree), alpha = 0.4) +
  geom_ribbon(aes(x = depth, ymin = ymin_tree, ymax = ymax_tree, fill = tree, group = tree), alpha = 0.3)+
  geom_line(aes(x = depth, y = mean_tree, col = tree, group = tree))+
  facet_wrap(~species, scales = "free", ncol = 4) +
  xlim(0,8) +
  theme_minimal() +
  labs(x = "Sensor depth (cm)", y = "Sap flux per section")

# ...individual level --------
ggplot(datap)+
  geom_point(aes(x = depth, y = flux, col = species), alpha = 0.4) + 
  geom_line(aes(x = depth, y = mean, col = species, group = interaction(tree, date)), alpha = 0.6)+
  geom_ribbon(aes(x = depth, ymin = ymin_tree, ymax = ymax_tree, fill = species), alpha = 0.3)+
  geom_line(aes(x = depth, y = mean_tree, col = species, group = tree))+
  facet_wrap(~paste(species, "-", tree), scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Sensor depth (cm)", y = "Sap flux per section")


# ...residual diagnostic plots -------
# residuals vs. fitted
datap %>% ggplot(aes(x = mean)) +
  geom_pointrange(aes(y = flux - mean,
                      ymin = flux - ymin,
                      ymax  =  flux - ymax),
                  alpha = 0.2) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_smooth(aes(y = flux - mean), se = F)
# some abberant measurement days produce outliers, but other than that
# the data are well covered

# predicted vs. observed values
datap %>% ggplot(aes(x = flux)) +
  geom_pointrange(aes(y = mean, ymin = ymin, ymax = ymax),
                  alpha = 0.2) +
  geom_abline(col = 2)  + 
  coord_flip() +
  theme_minimal() +
  labs(y = "Observed SFS", x = "Predicted SFS")
