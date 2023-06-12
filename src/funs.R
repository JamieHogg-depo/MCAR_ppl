# Functions

myfuns <- list()

# -----------------------------------------------------------------------------
#' @param W binary contiguity matrix (must be complete)
#' @param type (defaults to 'pcar') but also takes 'lcar'
myfuns$prep4MCAR <- function(W, type = "pcar"){
  if(type == "pcar"){
    
    # create sparse matrices
    W <- Matrix::Matrix(W, sparse = TRUE)
    I <- Matrix::Diagonal(nrow(W))
    D <- Matrix::Diagonal( x = Matrix::rowSums(W) )
    C <- solve(D) %*% W
    
    # Eigenvalues of C
    C_eigenvalues <- eigen(C)$values
    
    # get the CRS representation of C
    crs <- rstan::extract_sparse_parts(I+C)
    nC_w <- length(crs$w)
    
    # prepare output list
    return(
      list(C_eigenvalues = C_eigenvalues, 
           nC_w = nC_w,
           C_w = crs$w,
           C_v = crs$v,
           C_u = crs$u,
           D_id_C_w = which(crs$w == 1),
           offD_id_C_w = which(crs$w != 1))
    )
    
  }
  if(type == "lcar"){
    
    # create sparse matrices
    W <- Matrix::Matrix(W, sparse = TRUE)
    D <- Matrix::Diagonal( x = Matrix::rowSums(W) )
    I <- Matrix::Diagonal(nrow(W))
    C <- I - D + W
    # C and W only differ by the diagonal values
    # C has -1 on off diagonals
    
    # ISSUE: Diagonal element of C is zero if area has only one neighbor
    
    # get indices for diagonals
    jt <- rstan::extract_sparse_parts(W + 5*I) # 5 is arbritary
    # 5's will only be on the diagonals
    D_id_C_w <- which(jt$w == 5) 
    # any values that are not 5 are off diagonals
    offD_id_C_w <- which(jt$w == 1)
    
    # Eigenvalues of C
    C_eigenvalues <- eigen(C)$values
    
    # get the CRS representation of C
    # add an extra 1 to all diagonals to ensure they
    # are captured by `extract_sparse_parts`
    crs <- rstan::extract_sparse_parts(C + I)
    nC_w <- length(crs$w)
    
    # Remove 1 from the diagonals 
    crs$w[D_id_C_w] <- crs$w[D_id_C_w] - 1
    
    # prepare output list
    return(
      list(C = as.matrix(C),
           C_eigenvalues = C_eigenvalues, 
           nC_w = nC_w,
           C_w = crs$w,
           C_v = crs$v,
           C_u = crs$u,
           D_id_C_w = D_id_C_w,
           offD_id_C_w = offD_id_C_w)
    )
    
  }
  
}

# -----------------------------------------------------------------------------
#' @title fitLeroux_mc
#' @param formula fixed effect formula
#' @param W binary weight matrix
#' @param data data for the fixed effects
#' @param burnin number of iterations for burnin
#' @param n.sample number of iterations after burnin
#' @param n.chains number of chains; defaults to 4
#' @param family character; defaults to "poisson"
#' @param thin numeric; defaults to 1

myfuns$fitLeroux_mc <- function(formula, 
                         W, 
                         data, 
                         burnin, 
                         n.sample, 
                         n.chains = 4,
                         family = "poisson",
                         thin = 1){
  
  chains <- list()
  
  # fit 4 chains in series
  chains <- replicate(n.chains, 
                      S.CARleroux(as.formula(formula), 
                                  W = W, 
                                  data = data, 
                                  family = family, 
                                  thin = thin, 
                                  burnin = burnin, 
                                  n.sample = n.sample, 
                                  verbose = TRUE,
                                  prior.tau2 = c(1,0.5),
                                  prior.var.beta = 1), 
                      simplify = F)
  
  # get variable names
  var_names <- names(chains[[1]]$samples)
  # drop observed values from samples - assuming complete y
  var_names <- var_names[var_names != "Y"]
  
  # join all chains
  samples <- list()
  # loop through all variables
  for(i in 1:length(var_names)){
    # NOTE: all chains are the same model so we arbitrarily pick 
    # the first chain to extract variable names
    id <- which(names(chains[[1]]$samples) == var_names[i])
    #get dimensions for one chain
    cur_dims <- dim(chains[[1]]$samples[[id]])
    
    # collapse lists
    il <- list()
    for(j in 1:n.chains){
      # loop through each chain and add to il list
      il[[j]] <- chains[[j]]$samples[[id]]
    }
    samples[[var_names[i]]] <- do.call(rbind, il)
  }
  
  # add indicator for chains
  samples <- lapply(samples, function(x){
    data.frame(x) %>% 
      mutate(# integers from 1:n.chain
        .chain = sort(rep(1:n.chains, cur_dims[1])),
        # unique within each chain
        .iteration = rep(1:cur_dims[1], n.chains),
        # total number of iterations across all chains
        .draw = 1:(n.chains * cur_dims[1]) 
      )})
  
  # return samples
  return(samples)
  
}

## ----------------------------------------------------------------------------
#' @param samples output from fitLeroux_mc
#' @return dataset of summary of results including basic diagnostics

myfuns$summaryLeroux_mc <- function(samples){
  
  # get initial objects
  var_names <- names(samples)
  within_vars_length <- lengths(samples) - 3
  
  ll <- list()
  for(i in 1:length(var_names)){
    if(within_vars_length[i] > 1){
      ll[[i]] <- samples[[i]] %>% 
        posterior::summarise_draws() %>% 
        dplyr::mutate(variable = paste0(var_names[i], "[", 1:within_vars_length[i], "]"))
    }else{
      ll[[i]] <- samples[[i]] %>% 
        posterior::summarise_draws() %>% 
        dplyr::mutate(variable = var_names[i])
    }
  }
  
  # output summary dataset
  return(bind_rows(ll))
  
}