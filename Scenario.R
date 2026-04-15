library(DCL)
library(tidyverse)
library(magrittr)
library(emdist)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

source("~/Documents/DISSERTATION/Project/Data and Utils.R")
source("~/Documents/DISSERTATION/Project/Deterministic.R")
source("~/Documents/DISSERTATION/Project/Shock.R")
source("~/Documents/DISSERTATION/Project/Metrics.R")

#                   SCENARIO ASSIGNMENT

omega <- c(0, 0.5, 0.75, 1, 1.5)

y_list <- function(){
  #Creates a list of calendar years to introduce the shock
  #Generally half of m, then last 4 years
  #y = i + j - 1
  if(m < 4){
    y <- c(1:m)
  } else{
    y <- c(
      (m + 1) %/% 2,
      m-3,
      m-2,
      m-1,
      m
    )
  }
  sort(y)
}
y <- y_list()

scenario_matrix <- crossing(y = y, omega = omega) %>%
  filter(omega != 0 | y == min(y))

#                   RUN SCENARIO
run_scenario <- function(y, omega, N = N_det, X = X_det, I = I_det){
  #Brings together metrics and model fitting and runs for a single scenario
  shock_data <- apply_shock(X = X, I = I, y = y, omega = omega, N = N)
  N_true <- shock_data$N
  X_true <- shock_data$X
  I_true <- shock_data$I
  
  N_obs <- upper_triangle(N_true)
  X_obs <- upper_triangle(X_true)
  I_obs <- upper_triangle(I_true)
  
  fitted <- fit_models(N_obs, X_obs, I_obs)
  model_names <- names(fitted$pred)
  
  metrics_list <- vector("list", length(model_names))
  cy_list <- vector("list", length(model_names))
  names(cy_list) <- model_names
  rows <- vector("list", length(model_names))
  
  for (i in seq_along(model_names)){
    model_iter <- model_names[i]
    X_pred <- fitted$pred[[model_iter]]
    metrics_obj <- calc_metrics(X_true, X_pred)
    cy_err_metrics <- calc_cy_err(X_true, X_pred)
    rows[[i]] <- tibble(
      y = y,
      omega = omega,
      model = model_iter,
      metrics = list(metrics_obj),
      cy_rel = list(cy_err_metrics),
      fitted = list(fitted$obj[[model_iter]]),
      pred = list(fitted$pred[[model_iter]]),
      truth = list(list(N = N_true, X = X_true, I = I_true))
    )
  }
  bind_rows(rows)
}
