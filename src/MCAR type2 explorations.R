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

M <- 9 #3000
K <- 2
W <- getGRID(M)$W
#W_inv <- solve(W) # SINGULAR
M <- nrow(W)
n <- M*K

# Spatial precision matrix
I <- diag(x=1, nrow = M)
D <- diag(rowSums(W))
D_inv <- diag(1/rowSums(W))
rho <- c(0.8, 0.9)
rho_D <- diag(rho)

# Risk factor matrices
Omega_R <- matrix(c(1.2, 0.7,
                    0.7, 3), byrow = T, ncol = K)
Omega_R <- Omega_R %*% t(Omega_R)
Sum_R <- solve(Omega_R)
Sum_R_chol <- chol(Sum_R)

## Exploring P0349 -------------------------------------------------------------
# Get Q and triangle
D_inv_sqr = diag(rowSums(W)^(-1/2))
D_sqr = diag(rowSums(W)^(1/2))
temp <- eigen(D_inv_sqr %*% W %*% D_inv_sqr)
Q <- temp$vectors
Tri <- diag(temp$values)

# check decomposition
round(temp$vectors %*% diag(temp$values) %*% solve(temp$vectors), 3)
round(D_inv_sqr %*% W %*% D_inv_sqr, 3)

# derive quantities before 11
Om <- list(); Om_sqr <- list(); A <- list(); G_j <- list(); Tt <- list()
for(i in 1:2){
  Om[[i]] <- diag(1-rho[i]*temp$values)
  Om_sqr[[i]] <- diag(sqrt(1-rho[i]*temp$values))
  A[[i]] <- D_sqr %*% Q %*% Om_sqr[[i]] %*% t(Q)
  G_j[[i]] <- A[[i]] %*% solve(A[[i]])
  Tt[[i]] <- A[[i]] %*% t(A[[i]]) # equal to D - rho[i] * W
}
G <- Matrix::bdiag(G_j)

# Equation 11 
round(solve(G) %*% kronecker(Sum_R, Tt[[1]]) %*% t(solve(G)),2)
round(Omega_R[1,1] * Tt[[1]], 2)

## -----------------------------------------------------------------------------

# get A
A <- kronecker(D, Omega_R)
A_inv <- kronecker(D_inv, Sum_R)

# get B
B <- kronecker(W, Omega_R %*% rho_D %*% Omega_R)
B_inv <- kronecker(solve(W), solve(Omega_R %*% rho_D %*% Omega_R)) # SINGULAR

# get Omega_A
Omega_A <- A - B

# Woodbury matrix identity
Sum_A <- A_inv - A_inv %*% solve( A %*% -B_inv + diag(1, nrow = M*K) )

# determinant of Omega_A
det(Omega_A)

# let vec(ABC) = (C^T bigotimes A)vec(B)
# let B be a K times M matrix of all ones
as.numeric(Omega_R %*% rho_D %*% Omega_R %*% matrix(1, nrow = K, ncol = M) %*% t(W))
kronecker(W, Omega_R %*% rho_D %*% Omega_R) %*% as.numeric(rep(1,M*K))

## Matrix determinant lemma ## -------------------------------------------------

# Define two matrices A and B
svd_B <- svd(B)
U <- svd_B$u %*% diag(sqrt(svd_B$d))
V <- svd_B$v %*% diag(sqrt(svd_B$d))
temp <- determinant(diag(1, M*K) - t(V)%*%A_inv%*%U, logarithm = T)$modulus
#B; U %*% t(V)

# Compute the determinant of A-B using the formula
s = Sys.time(); (t = determinant(A - B, logarithm = T)$modulus); Sys.time()-s
#det(A - U %*% t(V))
s = Sys.time(); (t = temp + determinant(A, logarithm = T)$modulus); Sys.time()-s
s = Sys.time(); (t = temp + K*sum(log(diag(D))) + M * determinant(Omega_R, logarithm = T)$modulus); Sys.time()-s


s = Sys.time(); (t = temp*prod(diag(D))^K * det(Omega_R)^M); Sys.time()-s

# Singular value decomposition of B
svd_B <- svd(B)
round(svd_B$u, 4)

svd_W <- svd(W)

svd_2 <- svd(Omega_R %*% rho_D %*% Omega_R)

round(kronecker(svd_2$u, svd_W$u), 4)

round(svd_B$u %*% diag(sqrt(svd_B$d)), 4)
round(kronecker(svd_2$u %*% diag(sqrt(svd_2$d)), svd_W$u %*% diag(sqrt(svd_W$d))), 4)
          
## Exploring P0355 -------------------------------------------------------------

B <- matrix(c(0.8,0.4,
              0.4,0.1), nrow = K, ncol = K)
B <- matrix(c(0.99,0,
              0,0.99), nrow = K, ncol = K)
Ip <- diag(1, nrow = K)
round(eigen(kronecker(Ip, D) - kronecker(B, W))$values, 3)
Sum_U <- solve(kronecker(Ip, D) - kronecker(B, W))

# draw a random u
u <- rmvnorm(1, rep(0, M*K), Sum_U)

# convert to phi
A <- Sum_R_chol
kronecker(A, diag(1, nrow = M)) %*% matrix(u, nrow = K*M)
