# Fitting the MCAR in Stan

library(tidyverse)
library(sf)
library(spdep)
library(mvtnorm)
library(rstan)
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

#' @param input stan summary (summary(fit)$summary)
getSubsetSummary <- function(input, regex){
  as.data.frame(input) %>% 
    rownames_to_column("parameter") %>% 
    relocate(parameter) %>% 
    filter(str_detect(parameter, regex))
}

#' @param W binary contiguity matrix (must be complete)
prep4MCARt1 <- function(W){
  
  # create sparse matrices
  W <- Matrix::Matrix(W, sparse = TRUE)
  D <- Matrix::Diagonal( x = Matrix::rowSums(W) )
  I <- Matrix::Diagonal(nrow(W))
  
  C <- I - D + W
  
  # Eigenvalues of C
  C_eigenvalues <- eigen(C)$values
  
  # get the CRS representation of C
  crs <- rstan::extract_sparse_parts(C)
  nC_w <- length(crs$w)
  
  # prepare output list
  return(
    list(C_eigenvalues = C_eigenvalues, 
         nC_w = nC_w,
         C_w = crs$w,
         C_v = crs$v,
         C_u = crs$u,
         D_id_C_w = which(crs$w != 1),
         offD_id_C_w = which(crs$w == 1))
  )
}

## Generate data ## ------------------------------------------------------------

M <- 1000
K <- 3
W <- getGRID(M)$W
M <- nrow(W)
n <- M*K

# Spatial precision matrix
I <- diag(x=1, nrow = M)
D <- diag(rowSums(W))
C <- I - D + W
rho <- 0.9
Omega_S = I - rho * C
Omega_Ss <- as(Omega_S, "sparseMatrix")

# get bounds for rho
D_inv_sqr <- diag((rowSums(W))^(-1/2))
rho_bounds <- c(1/min(eigen(D_inv_sqr %*% W %*% D_inv_sqr)$values),
                1/max(eigen(D_inv_sqr %*% W %*% D_inv_sqr)$values))

# Eigenvalues
e <- eigen(C) # BOTTLENECK
Q_c <- e$vectors
C_eigenvalues <- e$values

# Risk factor matrices
Omega_R <- matrix(c(1.2, 0.7, 1,
                    0.7, 3,   1,
                    1,   1,   2), byrow = T, ncol = K)
Omega_R <- Omega_R %*% t(Omega_R)
Sigma_R <- solve(Omega_R)

Cor_R <- cov2cor(Sigma_R) # sampling this
Cor_R_chol <- t(chol(Cor_R)) # sampling this
Sigma_sd <- diag(sqrt(diag(Sigma_R)))

# Get Omega_R from Cor_R
diag(1/sqrt(diag(Sigma_R))) %*% solve(Cor_R) %*% diag(1/sqrt(diag(Sigma_R)))

# Get cholesky of Covariance matrix
Sigma_R_chol = t(chol(Sigma_R))
Sigma_R_chol_d = diag(Sigma_R_chol)

# FULL mvn matrix
Omega_A <- kronecker(Omega_S, Omega_R)
s=Sys.time()
Sigma_A <- solve(Omega_A) # BOTTLENECK: 3 mins for 6000x6000
Sys.time()-s

# Get random effects
set.seed(45)
y_v <- MASS::mvrnorm(1, rep(0,M*K), Sigma_A)
y_mat <- matrix(y_v, nrow = M, ncol = K, byrow = T)

## Compile models ## -----------------------------------------------------------

comp_n <- stan_model(file = "working_src/MCAR/MCAR_naive.stan")
comp_e <- stan_model(file = "working_src/MCAR/MCAR_eff.stan")
comp_ef <- stan_model(file = "working_src/MCAR/MCAR_stanform.stan")

## Naive ## --------------------------------------------------------------------

# data list
data <- list(M = M, K = K,
             y_vec = y_v, 
             y = y_mat, y_mat = y_mat, C_eigenvalues = C_eigenvalues,
             C = C, zero = rep(0, M*K), rho = 0.95)

# fit model
m_s <- Sys.time()
fit_n <- sampling(object = comp_n, 
                data = data, 
                chains = 4, iter = 2000, warmup = 1000,
                cores = 4, 
                pars = c("Omega_R", "Omega_A", "Omega_S"), 
                include = F)
(rtmins_n <- as.numeric(Sys.time() - m_s, units = "mins"))

# Summarise results
print(fit_n, digits_summary = 3)
summ_n <- summary(fit_n)$summary

## Efficient ## ----------------------------------------------------------------

# data list
data <- list(M = M, K = K,
             y_vec = y_v, 
             y = y_mat, y_mat = y_mat, C_eigenvalues = C_eigenvalues,
             C = C, zero = rep(0, M*K), rho = 0.95)

# fit model
m_s <- Sys.time()
fit_e <- sampling(object = comp_e, 
                  data = data, 
                  chains = 4, iter = 2000, warmup = 1000,
                  cores = 4, 
                  pars = c("Sigma_R", "Omega_R", "ldet_C", "Cor_R_chol", "A_R", "A_S", "Omega_A", "Omega_S"), 
                  include = F)
(rtmins_e <- as.numeric(Sys.time() - m_s, units = "mins"))

# Summarise results
print(fit_e, digits_summary = 3)
summ_e <- summary(fit_e)$summary

## Efficient - with function ## ------------------------------------------------

# data list
data <- list(M = M, K = K,
             y_mat = y_mat, C_eigenvalues = C_eigenvalues,
             C = C)

# fit model
m_s <- Sys.time()
fit_ef <- sampling(object = comp_ef, 
                  data = data, 
                  chains = 4, iter = 2000, warmup = 1000,
                  cores = 4, 
                  pars = c("Sigma_R", "Omega_R", "Omega_S"), 
                  include = F)
(rtmins_ef <- as.numeric(Sys.time() - m_s, units = "mins"))

# Summarise results
print(fit_ef, digits_summary = 3)
summ_ef <- summary(fit_ef)$summary

## Compare ## ------------------------------------------------------------------

# Identical inference
getSubsetSummary(summ_n, "ll_mcar|rho")
getSubsetSummary(summ_e, "ll_mcar|rho")
getSubsetSummary(summ_ef, "ll_mcar|rho")

# Run time 
rtmins_e
rtmins_n
rtmins_ef

# Comparison of sampling efficiency
# ESS per minute on average
summary(as.data.frame(summ_ef)$n_eff/rtmins_ef)
summary(as.data.frame(summ_e)$n_eff/rtmins_e)
summary(as.data.frame(summ_n)$n_eff/rtmins_n)

## Reconstruct matrices ## -----------------------------------------------------

temp <- summ %>% as.data.frame() %>% 
  dplyr::filter(str_detect(row.names(summ), "Omega_R")) %>% 
  dplyr::select(mean) %>% 
  rownames_to_column(var = "parameter") %>% 
  separate(parameter, into = c("name", "i", "j"), sep = "\\[|\\]|,") %>% 
  mutate(i = as.numeric(i),
         j = as.numeric(j))

Omega_R <- matrix(NA, nrow = 4, ncol = 4)
for(n in 1:nrow(temp)){
  Omega_R[temp$i[n], temp$j[n]] <- temp$mean[n]
}

solve(Omega_R)
print(fit, digits_summary = 3, pars = "Sigma_R")
