functions{
/**
* Log probability density of the multivariate proper conditional autoregressive (MCAR) model
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
    matrix y_mat,               // input pars in matrix form M x K
    real rho,                   // spatial smoothing parameter
    matrix Omega_R,             // Precision matrix for factors
    matrix Sigma_R,             // Covariance matrix for factors
    vector C_w , 
    int [] C_v , 
    int [] C_u , 
    int [] offD_id_C_w ,        // indexes for off diagonal terms
    int [] D_id_C_w ,       // indexes for diagonal terms - length M
    vector C_eigenvalues,       // eigenvalues for C
    vector D_W,
    int M,                      // Number of areas
    int K) {                    // Number of factors
        vector[M] ldet_C;
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
		ldet_C = log1m( rho * C_eigenvalues );
        //for(i in 1:M){
            // equivelant to log(1-rho * C_eigenvalues[i])
            //ldet_C[i] = log1m( rho * C_eigenvalues[i] ); 
        //}
    return -0.5 * ( 
        M*K*log( 2 * pi() ) 
        -  ( - 2 * M * sum( log_d_Sigma_R_chol ) + K * sum( ldet_C ) + K * sum( log(D_W)) ) + sq_term 
        );
}
}
data{
	int<lower=1> M; 				// number of areas
	int<lower=1> K; 				// dimension of health conditions
	vector[M] D_W;					// vector with number of neighbors
	matrix[M,K] theta_mat;
	// sparse components of C matrix
	vector[M] C_eigenvalues;
	int nC_w;
	vector[nC_w] C_w; 
	int C_v[nC_w]; 
	int C_u[M+1]; 
	int offD_id_C_w[nC_w - M];		// indexes for off diagonal terms
	int D_id_C_w[M]; 				// indexes for diagonal terms - length M
	real rho_set;
}
parameters{
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
  
  // Construct correlation matrix Lower cholesky
  // R = Lcorr * Lcorr'
	Cor_R = multiply_lower_tri_self_transpose(Cor_R_chol); 
  
  // quad_form_diag: diag_matrix(sig) * R * diag_matrix(sig)
	Sigma_R = quad_form_diag(Cor_R, sd_R); 
  
  // Get precision
	// equivalent to inverse(Sigma)
	Omega_R = quad_form_diag(chol2inv(Cor_R_chol), 1 ./ sd_R);
}
model{
	target += std_normal_lpdf(sd_R); // prior for sigma
	target += lkj_corr_cholesky_lpdf(Cor_R_chol | 2.0); // prior for cholesky factor of a correlation matrix
	target += uniform_lpdf(rho | 0, 1);
	if(rho_set == 1){
		target += MPCAR_lpdf(theta_mat | 1.0, Omega_R, Sigma_R, 
		C_w, C_v, C_u, offD_id_C_w, D_id_C_w, C_eigenvalues, D_W,
		M, K);
	}else if(rho_set == 0){
		target += MPCAR_lpdf(theta_mat | 0.0, Omega_R, Sigma_R, 
		C_w, C_v, C_u, offD_id_C_w, D_id_C_w, C_eigenvalues, D_W,
		M, K);
	}else{
		target += MPCAR_lpdf(theta_mat | rho, Omega_R, Sigma_R, 
		C_w, C_v, C_u, offD_id_C_w, D_id_C_w, C_eigenvalues, D_W,
		M, K);
	}
}



