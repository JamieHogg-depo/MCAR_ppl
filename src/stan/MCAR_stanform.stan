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
real MCARt1_lpdf(
	matrix x_mat,				// input pars in matrix form M x K
	real rho, 					// spatial smoothing parameter
	matrix Omega_R, 			// Precision matrix for factors
	matrix Sigma_R, 			// Covariance matrix for factors
	matrix Omega_S, 			// Precision matrix for areas (M x M)
	matrix C, 					// C matrix (M x M)
	vector C_eigenvalues,		// eigenvalues for C
	int M, 						// Number of areas
	int K) { 					// Number of factors
		vector[M] ldet_C;
		vector[K] log_d_Sigma_R_chol = log(diagonal(cholesky_decompose(Sigma_R)));
		matrix[K, M] A_R = Omega_R * x_mat';
		matrix[M, K] A_S = Omega_S * x_mat;
		real sq_term = sum( A_R .* A_S' );
		for(i in 1:M){
			// equivelant to log(1-rho * lambda[i])
			ldet_C[i] = log1m( rho * C_eigenvalues[i] ); 
		}
		return -0.5 * ( 
		M*K*log( 2 * pi() ) 
		+  ( 2 * M * sum( log_d_Sigma_R_chol ) - K * sum( ldet_C ) ) + sq_term 
		);
}
}
data {
  int<lower=1> M; // number of areas
  int<lower=1> K; // dimension of health conditions
  matrix[M,K] y_mat;
  vector[M] C_eigenvalues;
  matrix[M,M] C;
}
parameters {
  cholesky_factor_corr[K] Cor_R_chol; // Lower Cholesky of R
  vector<lower=0>[K] sd_R;
  real<lower=0,upper=1> rho;
}
transformed parameters {
// Naming System:
// *What*_*SorR*_*anyspecifis*
// Cor_R_chol: correlation matrix for risk factor matrix
// lower Cholesky matrix

  corr_matrix[K] Cor_R;
  cov_matrix[K] Sigma_R; 
  matrix[K, K] Omega_R;
  matrix[M, M] Omega_S;
  
  // Construct correlation matrix Lower cholesky
  // R = Lcorr * Lcorr'
  Cor_R = multiply_lower_tri_self_transpose(Cor_R_chol); 
  
  // quad_form_diag: diag_matrix(sig) * R * diag_matrix(sig)
  Sigma_R = quad_form_diag(Cor_R, sd_R); 
  
  // Get precision
	// equivalent to inverse(Sigma)
	Omega_R = quad_form_diag(chol2inv(Cor_R_chol), 1 ./ sd_R);
	// equivalent to I - rho *C
	Omega_S = add_diag( - rho * C, 1 );
}
model {
  sd_R ~ cauchy(0, 5); // prior for sigma
  Cor_R_chol ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix
  y_mat ~ MCARt1(rho, Omega_R, Sigma_R, Omega_S, C, C_eigenvalues, M, K);
}
generated quantities{
	real ll_mcar = MCARt1_lpdf( y_mat | rho, Omega_R, Sigma_R, Omega_S, C, C_eigenvalues, M, K );
}

