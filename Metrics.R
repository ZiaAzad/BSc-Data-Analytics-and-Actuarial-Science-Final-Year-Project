source("~/Documents/DISSERTATION/Project/Data and Utils.R")
source("~/Documents/DISSERTATION/Project/Deterministic.R")
source("~/Documents/DISSERTATION/Project/Shock.R")

library(DCL)
library(tidyverse)
library(magrittr)
library(emdist)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

#                   METRICS


#Accuracy
cell_error <- function(X_true, X_pred){
  #Calculates cell level error between the deterministic dataset and the
  #model calculated value as in Agbeko et al. (2015)
  X_t <- X_true[holdout_idx]
  X_p <- X_pred[holdout_idx]
  result <- sqrt( sum( (X_t - X_p) ^ 2, na.rm = TRUE) / sum(X_t ^ 2, na.rm = TRUE) )
  if (result <= tol){
    return(0)
  }
  result
}


cal_year_error <- function(X_true, X_pred){
  #Calendar Year Error as in Agbeko et al. (2015), similar but not equivalent
  #to a RMSE measure
  cy <- CY[holdout_idx]
  X_t <- X_true[holdout_idx]
  X_p <- X_pred[holdout_idx]
  X_T_cy <- tapply(X_t, cy, sum, na.rm = TRUE)
  X_P_cy <- tapply(X_p, cy, sum, na.rm = TRUE)
  
  sqrt( sum( (X_T_cy - X_P_cy) ^ 2) / sum(X_T_cy ^ 2))
} 




relative_total_error <- function(X_true, X_pred){
  #Bias in forecasted vs actual values. Same as total error but signed instead
  #Positive = Over-Reserved (More conservative)
  #Negative = Under-Reserved (More risky)
  R_t <- sum(X_true[holdout_idx], na.rm = TRUE)
  R_p <- sum(X_pred[holdout_idx], na.rm = TRUE)
  
  if (R_t <= tol){
    return(0)
  }
  (R_p - R_t) / R_t
} 


#Accuracy + Stability
cy_rel_err <- function(X_true, X_pred){
  #Calculates a vector of error in each CY
  cy <- CY[holdout_idx]
  X_t <- outstanding_by_cy(X_true)
  X_p <- outstanding_by_cy(X_pred)
  aligned <- align_cy(X_t, X_p)
  t <- aligned$a
  p <- aligned$b
  err_rel <- (p - t) / t
  err_rel[!is.finite(err_rel)] <- 0
  err_rel
}




worst_cumulative_shortfall <- function(X_true, X_pred){
  #Worst shortfall scaled by max of what the ultimate value is
  #max point at which model is short of expected cash flow
  #Normalised to be comparable to different cases
  #Negative Value implies under reserving and more liquidity risk
  X_t <- outstanding_by_cy(X_true)
  X_p <- outstanding_by_cy(X_pred)
  aligned <- align_cy(X_t, X_p)
  t <- cumsum(aligned$a)
  p <- cumsum(aligned$b)
  
  loss <- p - t
  min(loss) / max(t)
}



#Stability
worst_error <- function(X_true, X_pred){
  #Evaluates worst error
  err <- cy_rel_err(X_true, X_pred)
  worst_idx <- which.max( abs(err))
  c(
    WorstCYError = unname(err[worst_idx]),
    WorstCY = as.numeric(names(err)[worst_idx])
  )
} 



wasserstein_dist <- function(X_true, X_pred){
  #It penalises where errors occur in time series as well as the value
  #This is important as the size of the error is just as important as the timing
  #particularly in reserving where under-estimation can be a risk in any year
  prob_X_t <- convert_cf_dist(X_true)
  prob_X_p <- convert_cf_dist(X_pred)
  aligned <- align_cy(prob_X_t, prob_X_p)
  pmf_t <- aligned$a
  pmf_p <- aligned$b
  
  cy_labels <- as.numeric( names(pmf_t))
  dist_T <- cbind( as.numeric(pmf_t), cy_labels)
  dist_P <- cbind( as.numeric(pmf_p), cy_labels)
  emd(dist_T, dist_P, dist = "manhattan")
}

#                   METRIC CALCULATION AND STORAGE
calc_metrics <- function(X_true, X_pred){
  #Creates a list of each metric which is called in each scenario
  c(
    CellError = cell_error(X_true, X_pred), #Diagnostic
    CYErrorRMSE = cal_year_error(X_true, X_pred),
    RelTotalError = relative_total_error(X_true, X_pred),
    WorstShortfall = worst_cumulative_shortfall(X_true, X_pred),
    WorstCYError = worst_error(X_true, X_pred),
    WassersteinDistance = wasserstein_dist(X_true, X_pred)
  )
}

#EXTRA NOT USED
abs_total_error <- function(X_true, X_pred){
  #Unsigned total reserve error
  R_t <- sum(X_true[holdout_idx], na.rm = TRUE)
  R_p <- sum(X_pred[holdout_idx], na.rm = TRUE)
  if (R_t <= tol) return(0)
  abs(R_p - R_t) / R_t
}

sum_abs_cell_error <- function(X_true, X_pred){
  #Sum of absolute cell errors
  X_t <- X_true[holdout_idx]
  X_p <- X_pred[holdout_idx]
  sum(abs(X_t - X_p), na.rm = TRUE)
}

calc_cy_err <- function(X_true, X_pred){
  err <- cy_rel_err(X_true, X_pred)
  out <- cbind(
    CY = as.numeric(names(err)),
    value = as.numeric(err)
  )
  colnames(out) <- c("CY", "value")
  out
}

#Measures are split into 3 categories
#Category 1 --- Accuracy ---
#Cell Error -- CY Error -- Relative Total Error

#Category 2 --- Stability ---
#CY Relative Error -- Worst CY Error -- Parameter Plots

#Category 3 --- Timing ---
#Earth Mover's Distance -- CY Relative Error
