# MCAR explorations

library(tidyverse)
library(sf)
library(spdep)
library(mvtnorm)
rm(list = ls())

## Functions ## ----------------------------------------------------------------

#' @title getGRID
#' @description creates a sf object that is a grid of specified size
#' @param M target number of areas
#' @returns A list of two objects: sf is a data.frame/sf object with the geometry for the grid, W is the binary weight matrix based on queen contiguity

getGRID <- function(M){
  
  out <- list()
  
  # provides the height and length dimensions for the grid that are as close to the specified M as possible
  dims <- c(floor(sqrt(M)), floor(M/floor(sqrt(M))))
  
  df <- data.frame(lon = c(0, dims[1]), 
                   lat = c(0, dims[2])) %>% 
    st_as_sf(coords = c("lon", "lat")) %>% 
    st_bbox() %>% 
    st_as_sfc()
  
  out$sf <- sf::st_make_grid(df, square = T, cellsize = c(1, 1)) %>% 
    sf::st_sf() %>% 
    dplyr::mutate(id = row_number())
  
  out$W <- nb2mat(poly2nb(out$sf, queen = T), style="B") #binary
  message(paste0("Created an sf object with ", prod(dims), " rows"))
  return(out)
}

#' @title getSparseElements
#' @param A square matrix - should be sufficiently sparse
getSparseElements <- function(A){
  
  # Ensure matrix is sparse
  A <- as(A, "sparseMatrix")
  
  # cell values
  cell_value <- A@x # non zero elements - length: number of non-zeros
  
  # column ids
  c_id <- A@i+1 # integer of column indices - length: number of non-zeros
  ua <- A@p+1 # integer indices for cell_value/c_id which shows when a new row starts
  # ua has length: number of rows of A
  nu_zeros_per_row <- ua[-1] - ua[-length(ua)]
  
  # row ids
  r_id <- rep(1:nrow(A), nu_zeros_per_row)
  
  # indexes for diagonals
  cell_value_id_diagonals <- which(r_id == c_id)
  cell_value_id_nondiagonals <- which(r_id != c_id)
  
  # return the list
  out <- list(cell_value = cell_value,
              c_id = c_id, 
              r_id = r_id, 
              cell_value_id_diagonals = cell_value_id_diagonals,
              cell_value_id_nondiagonals = cell_value_id_nondiagonals)
  return(out)
}

## Generate data ## ------------------------------------------------------------

M <- 100 #3000
K <- 4
W <- getGRID(M)$W
M <- nrow(W)
n <- M*K

# Spatial precision matrix
I <- diag(x=1, nrow = M)
D <- diag(rowSums(W))
C <- I - D + W
rho <- 0.95
Omega_S = I - rho * C
Omega_Ss <- as(Omega_S, "sparseMatrix")
Omega_S_crs <- getSparseElements(Omega_S)
#Sum_S = solve(Omega_S)

# Eigenvalues
e <- eigen(C) # BOTTLENECK
Q_c <- e$vectors
lamba_C <- e$values

# Risk factor matrices
Omega_R <- matrix(c(1.2, 0.7, 1,   1,
                    0.7, 3,   1,   1.1,
                    1,   1,   2,   0.8,
                    1,   1.1, 0.8, 0.3), byrow = T, ncol = K)
Omega_R <- Omega_R %*% t(Omega_R)
Sum_R <- solve(Omega_R)

Cor_R <- cov2cor(Sum_R) # sampling this
Cor_R_chol <- t(chol(Cor_R)) # sampling this
Sum_sd <- diag(sqrt(diag(Sum_R)))

# Get Omega_R from Cor_R
diag(1/sqrt(diag(Sum_R))) %*% solve(Cor_R) %*% diag(1/sqrt(diag(Sum_R)))

# Get cholesky of Covariance matrix
Sum_R_chol = t(chol(Sum_R))
Sum_R_chol_d = diag(Sum_R_chol)

# FULL mvn matrix
Omega_A <- kronecker(Omega_S, Omega_R)
s=Sys.time()
Sum_A <- solve(Omega_A) # BOTTLENECK: 3 mins for 6000x6000
Sys.time()-s

# Get random effects
set.seed(45)
y_v <- rnorm(M*K)
y <- matrix(y_v, nrow = M*K, ncol = 1)
X <- matrix(y_v, nrow = K, ncol = M, byrow = FALSE)
y_mat <- t(X)

## Derive log-likelihood ## ----------------------------------------------------

# Using matrix normal distribution
cat("\014")

# standard log-likelihood
s=Sys.time(); dmvnorm(y_v, mean = rep(0, n), sigma = Sum_A, log = T); Sys.time()-s

# I1 - without fast determinant
s=Sys.time()
-(n/2)*log(2*pi) + 0.5*(determinant(Omega_A, logarithm = TRUE)$modulus) - 0.5 * (t(y) %*% kronecker(Omega_S, Omega_R) %*% y)
Sys.time()-s

# I2 - fast implementation
s=Sys.time()
-(n/2)*log(2*pi) + 0.5*(-2*M*sum(log(Sum_R_chol_d)) + K*sum(log(1-rho*lamba_C))) - 0.5 * (t(y) %*% kronecker(Omega_S, Omega_R) %*% y)
Sys.time()-s

# I3 - faster implementation
s=Sys.time()
-(n/2)*log(2*pi) + 0.5*(-2*M*sum(log(Sum_R_chol_d)) + K*sum(log(1-rho*lamba_C))) - 0.5 * (t(y) %*% kronecker(Omega_Ss, Omega_R) %*% y)
Sys.time()-s

# I4 - Use trace by borrowing ideas from Matrix Normal distribution
s=Sys.time()
A_S <- Omega_S %*% y_mat
A_R <- Omega_R %*% t(y_mat)
-(n/2)*log(2*pi) + 0.5*(-2*M*sum(log(Sum_R_chol_d)) + K*sum(log(1-rho*lamba_C))) - 0.5 * psych::tr(A_R %*% A_S)
Sys.time()-s

# I5 - use componentwise product
s=Sys.time()
A_S <- Omega_S %*% y_mat
A_R <- Omega_R %*% t(y_mat)
-0.5*( n*log(2*pi) + (2*M*sum(log(Sum_R_chol_d)) - K*sum(log(1-rho*lamba_C))) + sum(A_R * t(A_S)) )
Sys.time()-s

'Our method is about 2000 and 50 times faster when using Omega_Ss and Omega_S, respectively.
For M = 1980 and K = 3, 0.0055 secs (I5) secs, 0.013 (I3), 37 sec (I1).
For M = 3782 and K = 3, our sparse method takes 0.097 secs as opposed to 8.6 mins. 5300 times faster
For M = 3782 and K = 3, faster implementation is around 40 times faster than fast implementation
For M = 3782 and K = 3, even faster implementation is around 460 and 3.15 times faster than fast and faster implementation,
Even faster implementation is 42899 times faster than standard implementation

For M = 2652 and K = 4, 0.00375 secs (I5) secs, 0.022 (I3), 1.92 secs (I2), 3.98 mins (I1).
For M = 2970 and K = 4, 0.053 secs (I5) secs, 0.0287 (I3), 2.53 secs (I2), 5.71 mins (I1).
'

## Sparse matrix multiplication ## ---------------------------------------------

W <- as(getGRID(10)$W, "sparseMatrix")
cell_value <- W@x # non zero elements - length number of non-zeros
c_id <- W@i+1 # integer of column indices - length number of non-zeros
ua <- W@p+1 # integer indices for cell_value/c_id which shows when a new row starts
# ua has length equal to number of rows
nu_zeros_per_row <- ua[-1] - ua[-length(ua)]
r_id <- rep(1:9, nu_zeros_per_row)

# Use Omega_Ss
cell_value <- Omega_Ss@x # non zero elements - length number of non-zeros
c_id <- Omega_Ss@i+1 # integer of column indices - length number of non-zeros
ua <- Omega_Ss@p+1 # integer indices for cell_value/c_id which shows when a new row starts
# ua has length equal to number of rows
nu_zeros_per_row <- ua[-1] - ua[-length(ua)]
r_id <- rep(1:nrow(Omega_Ss), nu_zeros_per_row)
cell_value_id_diagonals <- which(r_id == c_id)
cell_value_id_nondiagonals <- which(r_id != c_id)

#' @param A Omega_S_crs
#' @param y numeric vector of length M*K
#' @param row_nu numeric for the row number of Omega_S
#' @param K number of risk factors
getDiagSum <- function(A, y, row_nu, K){
  
  start <- 1+ K*(row_nu - 1)
  end <- row_nu * K
  
  mlt <- A$cell_value[A$cell_value_id_diagonals][row_nu] 

  colSums(mlt * Omega_R * matrix(y[start:end], nrow = K, ncol = K, byrow = F))
}

#' @param A Omega_S_crs
#' @param y numeric vector of length M*K
#' @param row_nu numeric for the row number of Omega_S
#' @param K number of risk factors
getOffDiagSum <- function(A, y, row_nu, K){
  
  start <- 1+ K*(row_nu - 1)
  end <- row_nu * K
  
  mlt <- A$cell_value[A$cell_value_id_nondiagonals][1] 
  
  colSums(mlt * Omega_R * matrix(y[start:end], nrow = K, ncol = K, byrow = F))
}

# use lapply
s=Sys.time()
diag_elements <- lapply(1:M, function(x) getDiagSum(A = Omega_S_crs, y = y, row_nu = x, K = K))
offdiag_elements <- lapply(1:M, function(x) getOffDiagSum(A = Omega_S_crs, y = y, row_nu = x, K = K))

# derive t(y) %*% Omega_A by columns
out <- list()
for(c in 1:M){
  ll <- list()
  r_sel <- Omega_S_crs$r_id[which(Omega_S_crs$c_id == c)]
  for(r in r_sel){
    if(r == c){
      ll[[r]] <- diag_elements[[c]]
    }else{
      ll[[r]] <- offdiag_elements[[r]]
    }
  }
  out[[c]] <- Reduce("+", Filter(Negate(is.null), ll))
}
delta <- unlist(out); Sys.time() - s
s=Sys.time(); delta2 <- as.numeric((t(y) %*% kronecker(Omega_Ss, Omega_R))); Sys.time() - s
identical(round(delta, 6), round(delta2, 6))

## DEPREC ## -------------------------------------------------------------------

# Speed tests with last component of log-likelihood
s=Sys.time()
(t(y) %*% kronecker(Omega_S, Omega_R) %*% y) # v1
Sys.time()-s

# derive matrix normal values
p1 <- Omega_S %*% t(X)
p2 <- Omega_R %*% X

s=Sys.time()
psych::tr(p1 %*% p2) # v2
Sys.time()-s

s=Sys.time()
psych::tr(p2 %*% p1) # v3
Sys.time()-s # about 4.3 times faster than v4

s=Sys.time()
(t(y) %*% kronecker(Omega_Ss, Omega_R) %*% y) # v4
Sys.time()-s # about 2.46 times faster than v2

## DEPREC ## -------------------------------------------------------------------

# Working with block diagonal matrices
t(y) %*% Omega_A
colSums(rbind(-rho * Omega_R * matrix(c(y[3:4], y[3:4]), byrow = F, nrow = 2), 
              Omega_R * matrix(c(y[1:2], y[1:2]), byrow = F, nrow = 2)))

Omega_R

-rho * Omega_R * matrix(1, nrow = K, ncol = K) %*% diag(y[1:3])

# how to get indexes of the y vector based on row index of Omega_S
row_nu <- 2
(start <- 1+ K*(row_nu - 1))
(end <- row_nu * K)





