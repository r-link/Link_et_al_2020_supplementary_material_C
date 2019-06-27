################################################################################
#	
#              HFD data analysis with nlme   
#
#  Aim: generate and inspect preliminary nonlinear mixed models for sap 
#       flow profiles with R package nlme (used for HMC starting values)
#
# Author: Roman Link    Contact: roman.link@plant-ecology.de
#
################################################################################

# Sections:
# 1. Data loading and preparation
# 2. Nonlinear regression analysis with nls/nlme
# 3. Model checking and evaluation
# 4. Results of the full model
# 5. Export model objects

# load packages ----------------------------------------------------------------
# list of packages 
pkgs <-c("tidyverse",     # framework for data manipulation
         "nlme")          # nonlinear mixed models
         

# check for existence of packages and install if necessary
to_install<-pkgs[!(pkgs %in% installed.packages()[,1])]
if (length(to_install)>0)  for (i in seq(to_install)) install.packages(to_install[i])
# load all required packages 
for (i in seq(pkgs)) require(pkgs[i], character.only = T)


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
         ASI      = scale1(ASI)) 

# preparation for nlme-based models: scale flux by sd(flux) to avoid crashing models
data1 <- data %>% mutate(flux0 = flux,
                         flux = flux0 / sd(flux0)) 

# 2. Nonlinear regression analysis with nls/nlme ------------------------------
# define model equation and store as a formula object for convenience
# cf. Eqn. 4 in the main text and Eqn. A.1 in the digital supplement
# mu is retransformed from logit transformation, K and c are retransformed from log 
# transformation within the equation
beta1 <- flux ~ exp(log_c) *
                 dbeta(x = reldepth,  
                       shape1 = exp(logit_mu) / (1 + exp(logit_mu)) * exp(log_k),        # alpha = mu * K
                       shape2  = (1 - exp(logit_mu) / (1 + exp(logit_mu))) * exp(log_k)) # beta  = (1 - mu) * K
# for simplicity and to reduce numerical problems, the density function of the beta distribution is 
# called directly as dbeta() [function calls are allowed in the model equations in nls() and nlme()]

# define nlme control parameters (increased iteration numbers)
cont <- nlmeControl(maxIter = 9999,
                    niterEM = 9999,
                    pnlsMaxIter = 9999,
                    minScale = 0.000001,
                    msMaxIter = 9999,
                    msVerbose = TRUE)


# 2.0 - fit simple nls model to get starting values for nlme -----
mod0 <- nls(formula = beta1, 
            data = data1,
            start = list(log_c = - 1, logit_mu = 0.13, log_k = 0.05)) 
                     # starting values are just educated guesses

# inspect model summary
summary(mod0)


# 2.1 - fit model A -----
# the most simple model considered in the paper:
# no fixed effects, no tree-specific differences in curve shape
model_A <- nlme(model = beta1,
             data  = data1, 
             fixed = list(log_c ~ 1, 
                          logit_mu ~ 1, 
                          log_k    ~ 1),
             random = list(species =         logit_mu + log_k ~ 1,
                           tree    = log_c                    ~ 1,
                           date    = log_c                    ~ 1),
             start = list(fixed = c(coef(mod0))),
             control = cont)   

# 2.2 - fit model B -----
# no fixed effects, species and tree level RE for mu and k
# the model is fitted using update() to use the estimates of model_A 
# as starting values
model_B <- update(model_A,
                random = list(species =         logit_mu + log_k ~ 1,
                              tree    = log_c + logit_mu + log_k ~ 1,
                              date    = log_c                    ~ 1)
                )

# 2.3 - fit model C -----
# full random effects and and parameter regressions for mu and k 
model_C <- update(model_B,
                  fixed = list(log_c ~ 1, 
                               logit_mu ~ WD + height + ASI, 
                               log_k    ~ WD + height + ASI),
                  # fixed effects starting values have to be defined 
                  # (estimates of model B for the intercepts and 0 for the effects)
                  start = list(fixed = c(fixef(model_B)[1:2],0,0,0, 
                                         fixef(model_B)[3],0,0,0))
                  )


# 3. Model checking and evaluation --------------------------------------------
# residual diagnostics
# model A
plot(model_A)
qqnorm(model_A, abline = c(0, 1))
# very obvious residual patterning indicates the poor fit of the model

# model B
plot(model_B)
qqnorm(model_B, abline = c(0, 1))
# looks much better than model A, but strong patterns remain

# model C
plot(model_C)
qqnorm(model_C, abline = c(0, 1))
# looks better, but obvious patterning and long tails of the distribution of the 
# residuals (driven by badly explained measurement days) indicate suboptimal fit
# (the much better look of the residuals in the model fitted in STAN indicates
# that nlme might have got stuck at a local optimum)

# compare by LRT, AIC and BIC
anova(model_A, model_B, model_C) 
# full model is best in terms of AIC and logLik/LRT, but not BIC

# compare by pseudo-r-squared, root mean square error and mean absolute error
# (based on checkmod() function in the utility script)
checkmod(data1$flux, predict(model_A))
checkmod(data1$flux, predict(model_B))
checkmod(data1$flux, predict(model_C))
# r-squared roughly the same between model B and C, 
# RMSE and MAE marginally lower (!) for model without fixed effects

# compare observed vs. predicted values
plot(data1$flux ~ predict(model_A), col = scales::alpha("steelblue", .3), pch = 20)
points(data1$flux ~ predict(model_B), col = scales::alpha(2, .3), pch = 20)
points(data1$flux ~ predict(model_C), col = scales::alpha(1, .15), pch = 20)
abline(0, 1, lty = 2)
# the much better fit of the tree-level models is evident

# 4. Results of the full model ------------------------------------------------
# summary of full model
summary(model_C)
# model indicates significant height and growth effects on mu and significant
# WD and height effects on K
# however, large correlations in the random effects and the shrunk-to-zero species
# level random effects for mu indicate that there might be some trouble with the model
# - it could well be a false convergence 

# calculate conditional pseudo-R^2
var(avg1$pred3)/(var(avg1$pred3) + var(avg1$resid))
# calculate pseudo-R^2 on stem level
var(avg1$pred2)/(var(avg1$pred2) + var(avg1$flux - avg1$pred2))
# calculate spec-level pseudo-R^2
var(avg1$pred1)/(var(avg1$pred1) + var(avg1$flux - avg1$pred1))
# calculate marginal pseudo-R^2
var(avg1$pred0)/(var(avg1$pred0) + var(avg1$flux - avg1$pred0))

# visualize model predictions
data1 %>%
  # get predictions on different levels
  mutate(pred0 = predict(model_C, level = 0),
         pred1 = predict(model_C, level = 1),
         pred2 = predict(model_C, level = 2),
         pred3 = predict(model_C),
         resid = resid(model_C)) %>%
  # plot for each stem separately with ggplot2
  ggplot(aes(x = depth, y = flux0)) +
  geom_point(aes(col = species), alpha = 0.4) + 
  geom_line(aes(y = pred2 * sd(flux0), col = species)) +
  geom_line(aes(y = pred3 * sd(flux0), col = species, group = date), lty = 2) +
  facet_wrap(~paste(species, "-", tree), scales = "free") +
  labs(x = "Distance from cambium (cm)", y = "Sap flux per section") +
  theme_minimal() + 
  theme(legend.position = "bottom") 

# 5. Export model objects -----------------------------------------------------
#export output (for starting values for bayesian models)
save(list = c("model_A", "model_B", "model_C"), file = "output/nlme_model_beta_ABC.RData")
