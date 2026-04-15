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
source("~/Documents/DISSERTATION/Project/Scenario.R")

#This file contains the reference scenario analysis
#To be used as a diagnostic check of results

X_det
N_det
I_det

test_p <- dcl.estimation(X_upper, N_upper, Tables = FALSE, adj = 1)
testfit <- dcl.predict(test_p, N_upper, Model = 2, Tail = FALSE, Tables = FALSE)
testfit_bdcl_par <- bdcl.estimation(Xtriangle = X_upper, Ntriangle = N_upper, 
                                Itriangle = I_upper, Tables = FALSE, adj = 1)
testfit_bdcl <- dcl.predict(testfit_bdcl_par, N_upper, Model = 1, Tail = FALSE, Tables = FALSE)

X_out <- outstanding_by_cy(X_det)
testfit$Dtotal[1:18]
X_out

test_omega <- 0.5
cy <- row(X_det) + col(X_det) - 1
test_y <- 15
mult <- ifelse(CY < test_y, 0, test_omega)
mult_i <- ifelse(CY < (test_y ), 0, test_omega)
X_sk <- X_det * (1 + mult)
I_sk <- I_det * (1 + mult_i)
X_det
X_sk
obs_dt <- upper_triangle(X_det)
obs_sk <- upper_triangle(X_sk)
clm_base <- clm(obs_dt)
clm_sk <- clm(obs_sk)
clm_base
clm_sk

dcl_test_par <- dcl.estimation(upper_triangle(X_sk), Ntriangle = N_upper, Tables = FALSE, adj = 1)
dcl_test_fit <- dcl.predict(dcl_test_par, N_upper, Model = 2, Tail = FALSE, Tables = FALSE)
dcl_test_par
clm_test <- outstanding_by_cy(clm_sk$triangle.hat)
dcl_test_fit$Dtotal[1:18]
clm_test
outstanding_by_cy(dcl_test_fit$Xtotal)

dcl_test_par2 <- bdcl.estimation(upper_triangle(X_sk), Ntriangle = N_upper, Itriangle = upper_triangle(I_sk), Tables = FALSE, adj = 1)
bdcl_test_fit <- dcl.predict(dcl_test_par2, N_upper, Tail = FALSE, Tables = FALSE, Model = 2)
outstanding_by_cy(bdcl_test_fit$Xtotal)

outstanding_by_cy(X_sk)

cat("\nDETERMINISTIC TRUTH:\n")
print(outstanding_by_cy(X_det))
cat("\nDETERMINISTIC DCL:\n")
print(testfit$Dtotal[1:18])
cat("\nDETERMINISTIC BDCL:\n")
print(testfit_bdcl$Dtotal[1:18])

cat("\nSHOCKED TRUTH:\n")
print(outstanding_by_cy(X_sk))
cat("\nSHOCKED CLM:\n")
print(outstanding_by_cy(clm_sk$triangle.hat))
cat("\nSHOCKED DCL:\n")
print(dcl_test_fit$Dtotal[1:18])
cat("\nSHOCKED BDCL:\n")
print(bdcl_test_fit$Dtotal[1:18])

testfit_bdcl
dcl_test_par2

pars <- data.frame(
  dcl_inflat_ns = testfit_bdcl_par$inflat.DCL,
  bdcl_inflat_ns = testfit_bdcl_par$inflat,
  dcl_inflat_s = dcl_test_par2$inflat.DCL,
  bcl_inflat_s = dcl_test_par2$inflat
)
pars

metrics_clm_base <- calc_metrics(X_det, clm_base$triangle.hat) 
metrics_dcl_base <- calc_metrics(X_det, testfit$Xtotal) 
metrics_bdcl_base <- calc_metrics(X_det, testfit_bdcl$Xtotal) 

metrics_clm_sk <- calc_metrics(X_sk, clm_sk$triangle.hat)
metrics_dcl_sk <- calc_metrics(X_sk, dcl_test_fit$Xtotal)
metrics_bdcl_sk <- calc_metrics(X_sk, bdcl_test_fit$Xtotal)

floor_tol <- function(x) ifelse(abs(x) < tol, 0, x)
fmt <- function(x) formatC(floor_tol(as.numeric(x)), format = "f", digits = 4)

metrics_df <- data.frame(
  Metric = names(metrics_clm_base),
  CLM = fmt(metrics_clm_base),
  DCL = fmt(metrics_dcl_base),
  BDCL = fmt(metrics_bdcl_base)
)
metrics_df

metrics_df_sk <- data.frame(
  Metric = names(metrics_clm_base),
  CLM = fmt(metrics_clm_sk),
  DCL = fmt(metrics_dcl_sk),
  BDCL = fmt(metrics_bdcl_sk)
)
metrics_df_sk

testfit_bdcl$Dtotal[37]

