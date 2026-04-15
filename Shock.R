source("~/Documents/DISSERTATION/Project/Data and Utils.R")
source("~/Documents/DISSERTATION/Project/Deterministic.R")



#                   RAMP SHOCK FUNCTION
shock_triangle_ramp <- function(triangle, type, y, omega, L = 4){
  #Alternative example of a shock which increases in magnitude over 4 years
  yr <- if (type == "Incurred") (y - 2) else y
  pct_shock <- (CY - yr + 1) / L
  pct_shock <- pmin(pmax(0, pct_shock), 1) #Cap at 1
  
  tri_copy <- triangle
  tri_copy[!is.na(triangle)] <- triangle[!is.na(triangle)] * (1 + omega * pct_shock[!is.na(triangle)])
  tri_copy
}



shock_triangle <- function(triangle, type, y, omega){
  m <- nrow(triangle)
  cy <- row(triangle) + col(triangle) - 1
  
  multiplier <- ifelse(cy < y, 1, (1 + omega))
  out <- triangle * multiplier
  out
}



#                   APPLY SHOCK TO DATA
apply_shock <- function(X, I, y, omega, N){
  #This function applies a shock to the paid and incurred data
  #Counts is unchanged for a clearer comparison but is taken as a parameter
  #for easier data handling
  
 
  shocked <- (CY >= y) & !is.na(X)
  
  X_shock <- shock_triangle(X, "Paid", y, omega)
  I_shock <- shock_triangle(I, "Incurred", y, omega)
  list(N = N, X=X_shock, I=I_shock)
}
