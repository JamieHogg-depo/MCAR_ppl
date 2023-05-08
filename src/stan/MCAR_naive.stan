functions{
real standardMCAR_lpdf(
	vector y_vec,
	vector zero,
	matrix Omega_A){
	return multi_normal_prec_lpdf( y_vec | zero, Omega_A );
}
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
  matrix[K, K] Omega_R;
  matrix[M, M] Omega_S;
  matrix[M*K, M*K] Omega_A;
  
  // Construct correlation matrix Lower cholesky
  // R = Lcorr * Lcorr'
  Cor_R = multiply_lower_tri_self_transpose(Cor_R_chol);  
  
  // Get precision
	// equivalent to inverse(Sigma)
	Omega_R = quad_form_diag(chol2inv(Cor_R_chol), 1 ./ sd_R);
	// equivalent to I - rho *C
	Omega_S = add_diag( - rho * C, 1 );
  
  // Get MCAR precision matrix
	Omega_A = kronecker_prod(Omega_S, Omega_R);
}
model{
  sd_R ~ cauchy(0, 5); // prior for sigma
  Cor_R_chol ~ lkj_corr_cholesky(2.0); // prior for cholesky factor of a correlation matrix
  y_vec ~ multi_normal_prec( zero, Omega_A );
}
generated quantities{
	real ll_mcar = multi_normal_prec_lpdf( y_vec | zero, Omega_A );
}

