###############################################################################
#
#  utility functions to simplify scripts
#
###############################################################################

# handling of stan output  ----------------------------------------------------
# function to extract data from a subset of chains and return data.frame
extract_chain <- function(stan_fit, parnames, chains = 1:4) {
  chains_out <- extract(stan_fit, parnames, 
                        permute = FALSE) %>% as.array 
  
  if(length(chains) == 1) return(as.data.frame(chains_out[ , chains , ]))
  else return( chains_out[ , chains , ] %>% 
                 apply(2, as.data.frame) %>%
                 bind_rows )
}

# function for mean and credible intervals 
mci_stat <- function(x, range = 0.95) {
  x[is.infinite(x)] <- NA # remove infinite values (can happen in dweibull for low x)
  lower <- (1 - range) / 2
  upper <- 1 - lower
  c(ymin = as.numeric(quantile(x, lower, na.rm = T)),
    y = mean(x, na.rm = T),
    ymax = as.numeric(quantile(x, upper, na.rm = T)))
} 

# function for means and 95% credible intervals
mci <- function(x) {
  out <- c(quantile(x, 0.025), mean(x), quantile(x, 0.975))
  names(out) <- c("ymin", "y", "ymax")
  return(out)
}

# model evaluation  -----------------------------------------------------------
# kvalseth pseudo r-squared
rsq_kval <- function(y, yhat) 1 - var(y - yhat) / var(y)
# nakagawa style pseudo r-squared
rsq_naka <- function(y, yhat) var(yhat) / (var(yhat) + var(y - yhat))
# root mean square error
rmse <- function(y, yhat) sqrt(mean((y - yhat) ^2))
# mean absolute error
mae  <- function(y, yhat) mean(abs(y - yhat))
# predictive bias
bias <- function(y, yhat) mean((yhat - y))


checkmod <- function(y, yhat){
  return(c(rsq_kval = rsq_kval(y, yhat), rsq_naka = rsq_naka(y, yhat),
           rmse = rmse(y, yhat), mae = mae(y, yhat),
           cor = cor(y, yhat), bias = bias(y, yhat)))
}


# general ---------------------------------------------------------------------
# wrapper around scale function that drops all attributes
scale1 <- function(x) as.numeric(scale(x))
# in newer versions of the tidyverse, calling scale() inside tibbles can cause problems
# because of the changed coercion rules

# function for the expected value of radial sap flow profiles 
bfun <- function (x, mu, k){
  alpha <- mu * k
  beta  <- (1 - mu) * k
  out   <- dbeta(x, alpha, beta)
  out[x > 1] <- NA
  return(out)
}
