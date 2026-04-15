library(DCL)

data("NtriangleBDCL")
data("XtriangleBDCL")
data("ItriangleBDCL")
N_tri <- NtriangleBDCL
X_tri <- XtriangleBDCL
I_tri <- ItriangleBDCL

m <- nrow(N_tri)
CY <- row(N_tri) + col(N_tri) - 1
#Diagonals assigned with same integer - top left is current year
#For the code, first calendar year is 1

holdout_matrix <- CY > m #Boolean Matrix where True is the unobserved region
holdout_idx <- which(holdout_matrix, arr.ind = TRUE) #Indices of the triangle
phi <- 0
tol <- 1e-4

#                   UPPER TRIANGLE
upper_triangle <- function(triangle){
  #This is just to get the upper left triangle to use for testing
  #so that the functions cannot see the lower data and distort results
  m <- nrow(triangle)
  copy_tri <- triangle
  copy_tri[ row(copy_tri) + col(copy_tri) > (m+1) ] <- NA
  copy_tri
}




#                   ALIGN MATRICES FOR DISTRIBUTIONS
align_cy <- function(a, b){
  yrs <- union(names(a), names(b))
  ord <- order(as.numeric(yrs)) 
  a <- a[yrs]
  b <- b[yrs]
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  
  a <- a[ord]
  b <- b[ord]
  list(a = a, b = b)
}



#                   AGGREGATION BY CY
outstanding_by_cy <- function(triangle){
  cy <- CY[holdout_idx]
  holdout_region <- triangle[holdout_idx]
  out <- tapply(holdout_region, cy, sum, na.rm = TRUE)
  out <- out[order(as.numeric(names(out)))]
  out
}



#                   CUMULATIVE TO INCURRED
cum_to_inc <- function(cum_triangle) {
  inc_triangle <- matrix(0, nrow = m, ncol = m)
  
  for(i in 1:m){
    inc_triangle[i, 1] <- cum_triangle[i, 1]
        for(j in 2:m){
      inc_triangle[i, j] <- cum_triangle[i, j] - cum_triangle[i, j-1]
    }
  }
  
  return(inc_triangle)
}



#                   CF DISTRIBUTION TIMING
convert_cf_dist<- function(X){
  #Convert predictions into CFs and then into a distribution summing to 1
  cy <- CY[holdout_idx]
  X_cy <- outstanding_by_cy(X)
  prob <- X_cy / sum(X_cy, na.rm = TRUE)
  
  prob
}



#                   DCL PREDICT
forecast.dcl <- function(parameters, model, tail = FALSE, num.dec = 4){
  #DCL predictions using model 2
  #Produces m x m matrix
  #I used this as it was easier for me to understand than the package
  out <- matrix(0, nrow = m, ncol = m)
  
  alpha.N <- as.vector(parameters$alpha.N)
  beta.N <- as.vector(parameters$beta.N)
  pj <- as.vector(parameters$pj)
  mu <- as.vector(parameters$mu.adj)
  gamma.sev <- if (model == "BDCL") parameters$inflat else parameters$inflat.DCL
  sev <- as.vector(gamma.sev * mu)
  Nhat <- outer(alpha.N, beta.N, "*")
  
  for (i in 1:m){
    for (j in 1:m){
      delay <- 0
      
      for (l in 1:j){
        index <- j - l + 1
        delay <- delay + beta.N[index] * pj[l]
      }
      
      out[i, j] <- alpha.N[i] * sev[i] * delay
    }
  }
  round(out, num.dec)
}



#                   FITTING MODELS
fit_models <- function(N, X, I){
  fit_clm <- clm(X)
  fit_dcl <- bdcl.estimation(X, N, I, adj = 1, Tables = FALSE)
  #It uses bdcl but it is not used at all - just for the logic of the function to work
  fit_bdcl <- bdcl.estimation(X, N, I, adj = 1, Tables = FALSE)
  
  pred_dcl <- forecast.dcl(fit_dcl, model = "DCL")
  pred_bdcl <- forecast.dcl(fit_bdcl, model = "BDCL")
  
  list(
    pred = list(
      CLM = as.matrix(fit_clm$triangle.hat),
      DCL = as.matrix(pred_dcl),
      BDCL = as.matrix(pred_bdcl)
    ),
    obj = list(
      CLM = fit_clm,
      DCL = fit_dcl,
      BDCL = fit_bdcl
    )
  )
}
