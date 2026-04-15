library(DCL)

source("~/Documents/DISSERTATION/Project/Data and Utils.R")

params <- params_det <- bdcl.estimation(X_tri, N_tri, I_tri, adj = 1, Tables = FALSE)

N_det <- params_det$Nhat
N_upper <- upper_triangle(N_det)

X_det <- forecast.dcl(parameters = params_det, model = "DCL")
X_upper <- upper_triangle(X_det)

I_raw_cum <- get.cumulative(I_tri)
I_cum_clm <- clm(I_raw_cum)
I_cum_hat <- I_cum_clm$triangle.hat
I_det <- cum_to_inc(I_cum_hat)
I_upper <- upper_triangle(I_det) 
#This produced more reliable results than using clm() on the raw incremental data

m <- nrow(X_det)
CY <- row(X_det) + col(X_det) - 1
holdout_matrix <- CY > m
holdout_idx <- which(holdout_matrix, arr.ind = TRUE)