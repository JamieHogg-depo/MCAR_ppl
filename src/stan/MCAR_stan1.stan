data {
  int<lower=1> M; // number of areas
  int<lower=1> K; // dimension of health conditions
  array[M] vector[K] y; // observations: a list of M vectors (each has K elements)
  matrix[M,K] y_mat;
  vector[M] C_eigenvalues;
  matrix[M,M] Omega_S;
  vector[K] zero; 
}
parameters {
  cholesky_factor_corr[K] Cor_R_chol; // Lower Cholesky of R
  vector<lower=0>[K] sd_R;
}
transformed parameters {
// Naming System:
// *What*_*SorR*_*anyspecifis*
// Cor_R_chol: correlation matrix for risk factor matrix
// lower Cholesky matrix

  corr_matrix[K] Cor_R; // correlation matrix
  cov_matrix[K] Sigma_R; // VCV matrix
  matrix[K, K] Omega_R;
  
  // Construct correlation matrix Lower cholesky
  // R = Lcorr * Lcorr'
  Cor_R = multiply_lower_tri_self_transpose(Cor_R_chol); 
  
  // quad_form_diag: diag_matrix(sig) * R * diag_matrix(sig)
  Sigma_R = quad_form_diag(Cor_R, sd_R); 
  
  // Get precision
  // equivalent to inverse(Sigma)
  Omega_R = quad_form_diag(chol2inv(Cor_R_chol), 1 ./ sd_R);
}
model {
  sd_R ~ cauchy(0, 5); // prior for sigma
  Cor_R_chol ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix
  y ~ multi_normal( zero, Sigma_R );
}
generated quantities{
	vector[M] ldet_C;
	real fast;
	real slow;
	vector[K] log_d_Sigma_R_chol = log(diagonal(cholesky_decompose(Sigma_R)));
	matrix[K, M] A_R = Omega_R * y_mat';
	matrix[M, K] A_S = Omega_S * y_mat;
	real sq_term = sum( A_R .* A_S' ); 
	for(i in 1:M){
		ldet_C[i] = log1m( 0.95 * C_eigenvalues[i] ); // equivelant to log(1-rho * lambda[i])
	}
	fast = -0.5 * ( M*K*log( 2 * pi() ) +  ( 2 * M * sum( log_d_Sigma_R_chol ) - K * sum( ldet_C ) ) + sq_term );
	slow = multi_normal_lpdf( y | zero, Sigma_R );
}

