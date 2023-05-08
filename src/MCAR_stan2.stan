functions{
matrix kronecker_prod(matrix A, matrix B) {
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
    for (j in 1:n) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start = (j - 1) * q + 1;
      col_end = (j - 1) * q + q;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
}
}
data {
  int<lower=1> M; // number of areas
  int<lower=1> K; // dimension of health conditions
  array[M] vector[K] y; // observations: a list of M vectors (each has K elements)
  matrix[M,K] y_mat;
  vector[M*K] y_vec;
  vector[M] C_eigenvalues;
  matrix[M,M] C;
  vector[M*K] zero; 
  //real rho;
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

  corr_matrix[K] Cor_R; // correlation matrix
  cov_matrix[K] Sigma_R; // VCV matrix
  matrix[K, K] Omega_R;
  matrix[M, M] Omega_S;
  matrix[M*K, M*K] Omega_A;
  
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
  
  // Get MCAR precision matrix
	//Omega_A = kronecker_prod(Omega_S, Omega_R);
  
  // Log Likelihood of MCAR
	vector[M] ldet_C;
	real fast;
	vector[K] log_d_Sigma_R_chol = log(diagonal(cholesky_decompose(Sigma_R)));
	matrix[K, M] A_R = Omega_R * y_mat';
	matrix[M, K] A_S = Omega_S * y_mat;
	real sq_term = sum( A_R .* A_S' ); 
	for(i in 1:M){
		ldet_C[i] = log1m( rho * C_eigenvalues[i] ); // equivelant to log(1-rho * lambda[i])
	}
	fast = -0.5 * ( M*K*log( 2 * pi() ) +  ( 2 * M * sum( log_d_Sigma_R_chol ) - K * sum( ldet_C ) ) + sq_term );
}
model {
  sd_R ~ cauchy(0, 5); // prior for sigma
  Cor_R_chol ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix
  target += fast;
  //y_vec ~ multi_normal_prec( zero, Omega_A );
}
//generated quantities{
//	real slow = multi_normal_prec_lpdf( y_vec | zero, Omega_A );
//}

