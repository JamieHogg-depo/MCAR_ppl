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
  vector[M*K] theta_v;
  matrix[M,M] C;
  vector[M] D_W;
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
		Omega_S = diag_matrix(D_W) * add_diag( - rho * C, 1 );

	// Get MCAR precision matrix
		Omega_A = kronecker_prod(Omega_S, Omega_R);
}
model{
	target += std_normal_lpdf(sd_R); // prior for sigma
	target += lkj_corr_cholesky_lpdf(Cor_R_chol | 2.0); // prior for cholesky factor of a correlation matrix
	target += uniform_lpdf(rho | 0, 1);
	target += multi_normal_prec_lpdf( theta_v | zero, Omega_A );
}

