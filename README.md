# Multivariate CAR models in Stan

This project will explore efficient methods of fitting separable multivariate conditional autoregressive (MCAR) priors in Stan. Specifically, we introduce efficient implementations for the multivariate:

  - Proper CAR
  - Intrinsic CAR
  - Leroux CAR
  
We also show how each of these can be used in the univariate case.

## Multivariate Proper CAR: `MPCAR_lpdf`

```
functions{
/**
* Log probability density of the multivariate conditional autoregressive (MCAR) model
* @param y_mat input pars in matrix form (M x K)
* @param rho Spatial dependence parameter
* @param C matrix that represents spatial neighborhood 
* @param C_eigenvalues eigenvalues of C
* @param M number of areas
* @param K number of factors
**
@return Log probability density
*/
real MPCAR_lpdf(
	matrix y_mat,				    // input pars in matrix form M x K
	real rho, 					    // spatial smoothing parameter
	matrix Omega_R, 			  // Precision matrix for factors
	matrix Sigma_R, 			  // Covariance matrix for factors
	vector C_w , 
	int [] C_v , 
	int [] C_u , 
	int [] offD_id_C_w ,		// indexes for off diagonal terms
	int [] D_id_C_w , 		  // indexes for diagonal terms - length M
	vector C_eigenvalues,		// eigenvalues for C
	vector D_W,
	int M, 						      // Number of areas
	int K) { 					      // Number of factors
		vector[M] diag_term;
		vector[K] log_d_Sigma_R_chol = log(diagonal(cholesky_decompose(Sigma_R)));
		matrix[K, M] A_R = Omega_R * y_mat';
		// Alternative specification for A_S
		// Omega_S as sparse matrix multiplied by columns of y_mat 
		// using crs_matrix_times_vector
			vector [ num_elements(C_w) ] ImrhoC;
			matrix[M, K] A_S;
			// Multiple off-diagonal elements by rho
			ImrhoC [ offD_id_C_w ] = - rho * C_w[ offD_id_C_w ];
			// Calculate diagonal elements of ImrhoC
			ImrhoC [ D_id_C_w ] = C_w[ D_id_C_w ];
			for(k in 1:K){
				A_S[,k] = csr_matrix_times_vector( M, M, ImrhoC, C_v, C_u, y_mat[,k] );
			}
		for(i in 1:M){
			diag_term[i] = D_W[i] * dot_product( A_S[i,], A_R[,i]);
		}  
    real sq_term = sum( diag_term );
    vector[M] ldet_C = log1m( rho * C_eigenvalues );
    return -0.5 * ( 
		M*K*log( 2 * pi() ) 
		-  ( - 2 * M * sum( log_d_Sigma_R_chol ) + K * sum( ldet_C ) + K * sum( log(D_W)) ) + sq_term 
		);
}
}
```

## Multivariate Intrinsic CAR: `MICAR_lpdf`

```
functions{
/**
* Log probability density of the multivariate conditional autoregressive (MCAR) model - type1
* @param y_mat input pars in matrix form (M x K)
* @param rho Spatial dependence parameter
* @param C matrix that represents spatial neighborhood 
* @param C_eigenvalues eigenvalues of C
* @param M number of areas
* @param K number of factors
**
@return Log probability density
*/
real MICAR_lpdf(
	matrix y_mat,				    // input pars in matrix form M x K
	matrix Omega_R, 			  // Precision matrix for factors
	vector C_w , 
	int [] C_v , 
	int [] C_u , 
	int [] offD_id_C_w ,		// indexes for off diagonal terms
	int [] D_id_C_w , 		  // indexes for diagonal terms - length M
	vector D_W,
	int M, 						      // Number of areas
	int K) { 					      // Number of factors
		vector[M] diag_term;
		matrix[K, M] A_R = Omega_R * y_mat';
		// Alternative specification for A_S
		// Omega_S as sparse matrix multiplied by columns of y_mat 
		// using crs_matrix_times_vector
			vector [ num_elements(C_w) ] ImrhoC;
			matrix[M, K] A_S;
			// Multiple off-diagonal elements by rho
			ImrhoC [ offD_id_C_w ] = - C_w[ offD_id_C_w ];
			// Calculate diagonal elements of ImrhoC
			ImrhoC [ D_id_C_w ] = C_w[ D_id_C_w ];
			for(k in 1:K){
				A_S[,k] = csr_matrix_times_vector( M, M, ImrhoC, C_v, C_u, y_mat[,k] );
			}
		for(i in 1:M){
			diag_term[i] = D_W[i] * dot_product( A_S[i,], A_R[,i]);
		}  
    real sq_term = sum( diag_term );
    return -0.5 * ( 
		  M*K*log( 2 * pi() ) + sq_term 
		);
}
}
```

## Multivariate Leroux CAR: `MLCAR_lpdf`

```
functions{
/**
* Log probability density of the multivariate conditional autoregressive (MCAR) model - type1
* @param x_mat input pars in matrix form (M x K)
* @param mu Mean vector
* @param tau Scale parameter - sigma^(-2)
* @param rho Spatial dependence parameter
* @param C matrix that represents spatial neighborhood 
* @param C_eigenvalues eigenvalues of C
* @param M number of areas
* @param K number of factors
**
@return Log probability density
*/
real MLCAR_lpdf(
	matrix x_mat,				// input pars in matrix form M x K
	real rho, 					// spatial smoothing parameter
	matrix Omega_R, 			// Precision matrix for factors
	matrix Sigma_R, 			// Covariance matrix for factors
	vector C_w , 
	int [] C_v , 
	int [] C_u , 
	int [] offD_id_C_w ,		// indexes for off diagonal terms
	int [] D_id_C_w , 		// indexes for diagonal terms - length M
	vector C_eigenvalues,		// eigenvalues for C
	int M, 						// Number of areas
	int K) { 					// Number of factors
		vector[M] ldet_C;
		vector[K] log_d_Sigma_R_chol = log(diagonal(cholesky_decompose(Sigma_R)));
		matrix[K, M] A_R = Omega_R * x_mat';
		// Alternative specification for A_S
		// Omega_S as sparse matrix multiplied by columns of x_mat 
		// using crs_matrix_times_vector
			vector [ num_elements(C_w) ] ImrhoC;
			matrix[M, K] A_S;
			// Multiple off-diagonal elements by rho
			ImrhoC [ offD_id_C_w ] = - rho * C_w[ offD_id_C_w ];
			// Calculate diagonal elements of ImrhoC
			ImrhoC [ D_id_C_w ] = 1 - rho * C_w[ D_id_C_w ];
			for(k in 1:K){
				A_S[,k] = csr_matrix_times_vector( M, M, ImrhoC, C_v, C_u, x_mat[,k] );
			}
		ldet_C = log1m( rho * C_eigenvalues );
		real sq_term = sum( A_R .* A_S' );
		return -0.5 * ( 
		M*K*log( 2 * pi() ) 
		+  ( 2 * M * sum( log_d_Sigma_R_chol ) - K * sum( ldet_C ) ) + sq_term 
		);
}
}
```
