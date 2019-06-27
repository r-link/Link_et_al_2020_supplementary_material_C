################################################################################
#	
#             Automated computation of starting values   
#
# Aim: generate starting values for non-linear hierarchical for vertical 
#      sap flow profiles (modified and randomly perturbed, starting from 
#      maximum likelihood estimate based on R package nlme)
#
# Author: Roman Link    Contact: rlink@gwdg.de
#
################################################################################

# The first part of this script (1-3) takes the model objects generated in the script
# "1_HFD_analysis_nlme.R", extracts all relevant model components and performs
# the transformations necessary to use them as starting values for the models in
# the "stan/" folder

# The second part of this script (4a-c) defines functions that return randomized 
# perturbations of the model components extracted before for the use as starting values
# in rstan::stan()

# WARNING: This script is only meant to be called internally from the script 
#          "2_HFD_analysis_stan.R" 
#          The required packages and dataset are loaded in this script

# -----------------------------------------------------------------------------

# Sections:

# 1. Starting values for model A
# 2. Starting values for model B
# 3. Starting values for model C
# 4. Functions for starting values 

# get final nlme models
load("output/nlme_model_beta_ABC.RData")

# 1. Starting values for model A --------------------------------------------
summ <- summary(model_A)

# a) starting values for fixed effects --------------------------------------
mult0_sc <- fixef(model_A)[1]
mu0  <- fixef(model_A)[2]
k0 <- fixef(model_A)[3]

# b) starting values for variance components---------------------------------
# raw variance
sigma_sc <- summ$sigma 

# c) starting values for hierarchical variance parameters-------------------

# c1) RE correlation matrices
# multivariate correlations between effects are modeled by decomposition
# of covariance matrix in scale vector and correlation matrix
speccormat <- corMatrix(model_A$modelStruct[[1]])$species
# the correlation matrices are cholesky-decomposed to improve model performance
# cholesky-decomposed correlation matrices
L_Omega_spec <- chol(speccormat)

# c2) RE standard deviations
# # day level standard deviation
tau_u <- as.numeric(VarCorr(model_A)[7,2])
# tree level standard deviations (reordered to correct sequence)
tau_v <- as.numeric(VarCorr(model_A)[c(5),2])
# species level standard deviations
tau_w <- as.numeric(VarCorr(model_A)[2:3,2])

# extract species level random effects matrix
spec_RE_matrix <- as.data.frame(ranef(model_A)[[1]]) %>%
  as.matrix

# calculate transposed random effects normalized with their standard deviations
# (starting values for scaled random effects are needed because raw REs are 
# pre-multiplied with their standard deviations in the model)
z_spec <- t(spec_RE_matrix %*% diag(1 / tau_w))

# get random effects for observation dates
u <- ranef(model_A)[[3]][,1] 
v <- ranef(model_A)[[2]][,1] 

# scale approximate REs
u_raw <- u/sd(u)
v_raw <- v/sd(v)

# get starting value for hierarchical scale parameter
hier_scale <- 1

# create list of starting values for models with QR & cholesky decomposition
start_A <- list(mu0 = mu0, 
                k0  = k0,
                mult0_sc = mult0_sc, 
                sigma_sc = sigma_sc, 
                u_raw = u_raw, 
                v_raw = v_raw, 
                z_spec = z_spec, 
                L_Omega_spec = L_Omega_spec, 
                tau_u = tau_u, 
                tau_v = tau_v, 
                tau_w = tau_w, 
                hier_scale = hier_scale, 
                spec_RE_matrix = spec_RE_matrix)


# 2. Starting values for model B --------------------------------------------
summ <- summary(model_B)

# a) starting values for fixed effects --------------------------------------
mult0_sc <- fixef(model_B)[1]
mu0  <- fixef(model_B)[2]
k0 <- fixef(model_B)[3]

# b) starting values for variance components---------------------------------
# raw variance
sigma_sc <- summ$sigma 

# c) starting values for hierarchical variance parameters-------------------

# c1) RE correlation matrices
# multivariate correlations between effects are modeled by decomposition
# of covariance matrix in scale vector and correlation matrix
treecormat <- corMatrix(model_B$modelStruct[[1]])$tree[c(2,3,1), c(2,3,1)]
speccormat <- corMatrix(model_B$modelStruct[[1]])$species
# the correlation matrices are cholesky-decomposed to improve model performance
# cholesky-decomposed correlation matrices
L_Omega_tree <- chol(treecormat)   
L_Omega_spec <- chol(speccormat)

# c2) RE standard deviations
# # day level standard deviation
tau_u <- as.numeric(VarCorr(model_B)[9,2])
# tree level standard deviations (reordered to correct sequence)
tau_v <- as.numeric(VarCorr(model_B)[c(6,7,5),2])
# species level standard deviations
tau_w <- as.numeric(VarCorr(model_B)[2:3,2])

# extract tree level random effects matrix
tree_RE_matrix <-
  as.data.frame(ranef(model_B)[[2]]) %>% 
  mutate(xx = rownames(.)) %>%
  separate(xx, into = c("spec", "tree"), sep = "/") %>%
  arrange(tree) %>% 
  select(-spec, -tree) %>%
  select(2,3,1) %>%
  as.matrix

# extract species level random effects matrix
spec_RE_matrix <- as.data.frame(ranef(model_B)[[1]]) %>%
  as.matrix

# calculate transposed random effects normalized with their standard deviations
# (starting values for scaled random effects are needed because raw REs are 
# pre-multiplied with their standard deviations in the model)
z_tree <- t(tree_RE_matrix %*% diag(1 / tau_v))
z_spec <- t(spec_RE_matrix %*% diag(1 / tau_w))

# get random effects for observation dates
u <- ranef(model_B)[[3]][,1] 

# scale approximate REs
u_raw <- u/sd(u)

# get starting value for hierarchical scale parameter
hier_scale <- 1

# create list of starting values for models with QR & cholesky decomposition
start_B  <- list(mu0 = mu0, 
                 k0  = k0,
                 mult0_sc = mult0_sc, 
                 sigma_sc = sigma_sc, 
                 u_raw = u_raw, 
                 z_tree = z_tree, 
                 z_spec = z_spec, 
                 L_Omega_tree = L_Omega_tree, 
                 L_Omega_spec = L_Omega_spec, 
                 tau_u = tau_u, 
                 tau_v = tau_v, 
                 tau_w = tau_w, 
                 hier_scale = hier_scale, 
                 tree_RE_matrix = tree_RE_matrix,
                 spec_RE_matrix = spec_RE_matrix)

# 3. Starting values for model C --------------------------------------------
summ <- summary(model_C)

# a) starting values for fixed effects --------------------------------------
# # untransformed parameter vectors
beta_mu0  <- fixef(model_C)[2:5]
beta_k0 <- fixef(model_C)[6:9]

# the transformed regression parameter vectors theta are calculated from beta_ß 
# as the corresponding parameter vector for the QR decomposition of the model matrix
# components of the QR decomposition
QR <- qr(data_fit$X)                   # qr decomposition of the model matrix
R_ast <- QR$qr[1:4,] * upper.tri(diag(4), diag = T) # corresponding R* matrix
# 
# # parameter vector for the regression for the shape parameter
theta_mu <- - as.numeric((R_ast / sqrt(data_fit$I - 1)) %*% beta_mu0)
# # parameter vector for the regression for the scale parameter
theta_k <- - as.numeric((R_ast / sqrt(data_fit$I - 1)) %*% beta_k0)

# average of multiplicator parameter
mult0_sc <- fixef(model_C)[1]

# b) starting values for variance components---------------------------------
# raw variance
sigma_sc <- summ$sigma 

# c) starting values for hierarchical variance parameters-------------------
# c1) RE correlation matrices
# multivariate correlations between effects are modeled by decomposition
# of covariance matrix in scale vector and correlation matrix
treecormat <- corMatrix(model_C$modelStruct[[1]])$tree[c(2,3,1), c(2,3,1)]
speccormat <- corMatrix(model_C$modelStruct[[1]])$species
# the correlation matrices are cholesky-decomposed to improve model performance
# cholesky-decomposed correlation matrices
L_Omega_tree <- chol(treecormat)   
L_Omega_spec <- chol(speccormat)

# c2) RE standard deviations
# # day level standard deviation
tau_u <- as.numeric(VarCorr(model_C)[9,2])
# tree level standard deviations (reordered to correct sequence)
tau_v <- as.numeric(VarCorr(model_C)[c(6,7,5),2])
# species level standard deviations
tau_w <- as.numeric(VarCorr(model_C)[2:3,2])

# extract tree level random effects matrix
tree_RE_matrix <-
  as.data.frame(ranef(model_C)[[2]]) %>% 
  mutate(xx = rownames(.)) %>%
  separate(xx, into = c("spec", "tree"), sep = "/") %>%
  arrange(tree) %>% 
  select(-spec, -tree) %>%
  select(2,3,1) %>%
  as.matrix

# extract species level random effects matrix
spec_RE_matrix <- as.data.frame(ranef(model_C)[[1]]) %>%
  as.matrix

# calculate transposed random effects normalized with their standard deviations
# (starting values for scaled random effects are needed because raw REs are 
# pre-multiplied with their standard deviations in the model)
z_tree <- t(tree_RE_matrix %*% diag(1 / tau_v))
z_spec <- t(spec_RE_matrix %*% diag(1 / tau_w))

# get random effects for observation dates
u <- ranef(model_C)[[3]][,1] 

# scale approximate REs
u_raw <- u/sd(u)

# get starting value for hierarchical scale parameter
hier_scale <- 1

# create list of starting values for models with QR & cholesky decomposition
start_B  <- list(theta_mu = theta_mu, 
                 theta_k  = theta_k,
                 mult0_sc = mult0_sc, 
                 sigma_sc = sigma_sc, 
                 u_raw = u_raw, 
                 z_tree = z_tree, 
                 z_spec = z_spec, 
                 L_Omega_tree = L_Omega_tree, 
                 L_Omega_spec = L_Omega_spec, 
                 tau_u = tau_u, 
                 tau_v = tau_v, 
                 tau_w = tau_w, 
                 hier_scale = hier_scale, 
                 tree_RE_matrix = tree_RE_matrix,
                 spec_RE_matrix = spec_RE_matrix)

# 4. Functions for starting values --------------------------------------------
# ---a) function that randomizes regression and variance parameters -----
#       (+- 10 percent of the estimated value)
randfun <- function(x) runif(length(x), 
                             min = min(c(0.9 * x, 1.1 * x)),
                             max = max(c(0.9 * x, 1.1 * x)))


# ---b) function that randomizes random effects -----
#       (raw random effects are scaled with SD = 1, 
#       so 0.1 is 10% of the random effects SD)
REfun   <- function(x) {
  if(is.matrix(x)){ 
    xval <- as.numeric(x)
    return(matrix(rnorm(length(xval), xval, 0.1), nrow = nrow(x)))
  }
  else return(rnorm(length(x), x, 0.1))
}


# ---c) functions that generate lists of randomized starting values ----
# because of problems when directly randomizing correlation matrices 
# (some turned out non positive definite, thus creating impossible starting
# values that caused Stan to crash) the REs are randomized at the untransformed
# scale and then back-transformed

# function for model A
startfun_A <- function(){
  # get randomized REs before cholesky decomposition 
  # (to get randomized positive definite correlation matrices)
  tau_w1 <- randfun(start_A$tau_w)
  spec_RE_matrix1 <- REfun(start_A$spec_RE_matrix)
  
  # prepare output
  out <- with(start_A, list(
    mu0  = randfun(k0),
    k0 = randfun(k0),
    mult0_sc = randfun(mult0_sc),
    sigma_sc = randfun(sigma_sc),
    u_raw = REfun(u_raw),
    v_raw = REfun(v_raw),
    z_spec = t(spec_RE_matrix1 %*% diag(1 / tau_w1)),
    L_Omega_spec = chol(cor(spec_RE_matrix1)),
    tau_u = randfun(tau_u),  
    tau_v = randfun(tau_v),
    tau_w = tau_w1,
    hier_scale = randfun(hier_scale)
  ))
  return(out)
}

# function for model B
startfun_B <- function(){
  # get randomized REs before cholesky decomposition 
  # (to get randomized positive definite correlation matrices)
  tau_v1 <- randfun(start_B$tau_v)
  tau_w1 <- randfun(start_B$tau_w)
  tree_RE_matrix1 <- REfun(start_B$tree_RE_matrix)
  spec_RE_matrix1 <- REfun(start_B$spec_RE_matrix)
  
  # prepare output
  out <- with(start_B, list(
    mu0  = randfun(k0),
    k0 = randfun(k0),
    mult0_sc = randfun(mult0_sc),
    sigma_sc = randfun(sigma_sc),
    u_raw = REfun(u_raw),
    z_tree = t(tree_RE_matrix1 %*% diag(1 / tau_v1)),
    z_spec = t(spec_RE_matrix1 %*% diag(1 / tau_w1)),
    L_Omega_tree = chol(cor(tree_RE_matrix1)),
    L_Omega_spec = chol(cor(spec_RE_matrix1)),
    tau_u = randfun(tau_u),  
    tau_v = tau_v1,
    tau_w = tau_w1,
    hier_scale = randfun(hier_scale)
  ))
  return(out)
}

# function for model C
startfun_C <- function(){
  # get randomized REs before cholesky decomposition 
  # (to get randomized positive definite correlation matrices)
  tau_v1 <- randfun(start_C$tau_v)
  tau_w1 <- randfun(start_C$tau_w)
  tree_RE_matrix1 <- REfun(start_C$tree_RE_matrix)
  spec_RE_matrix1 <- REfun(start_C$spec_RE_matrix)
  
  # prepare output
  out <- with(start_C, list(
    theta_mu  = randfun(theta_mu),
    theta_k = randfun(theta_k),
    mult0_sc = randfun(mult0_sc),
    sigma_sc = randfun(sigma_sc),
    u_raw = REfun(u_raw),
    z_tree = t(tree_RE_matrix1 %*% diag(1 / tau_v1)),
    z_spec = t(spec_RE_matrix1 %*% diag(1 / tau_w1)),
    L_Omega_tree = chol(cor(tree_RE_matrix1)),
    L_Omega_spec = chol(cor(spec_RE_matrix1)),
    tau_u = randfun(tau_u),  
    tau_v = tau_v1,
    tau_w = tau_w1,
    hier_scale = randfun(hier_scale)
  ))
  return(out)
}



