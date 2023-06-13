functions{
/**
* Log probability density of the multivariate leroux CAR model
* @param x_mat input pars in matrix form (M x K)
* @param rho Spatial dependence parameter
* @param C matrix that represents spatial neighborhood 
* @param C_eigenvalues eigenvalues of C
* @param M number of areas
* @param K number of factors
**
@return Log probability density
*/
real MLCAR_lpdf(
    matrix x_mat,               // input pars in matrix form M x K
    real rho,                   // spatial smoothing parameter
    matrix Omega_R,             // Precision matrix for factors
    matrix Sigma_R,             // Covariance matrix for factors
    vector C_w , 
    int [] C_v , 
    int [] C_u , 
    int [] offD_id_C_w ,        // indexes for off diagonal terms
    int [] D_id_C_w ,       // indexes for diagonal terms - length M
    vector C_eigenvalues,       // eigenvalues for C
    int M,                      // Number of areas
    int K) {                    // Number of factors
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

data{
	int<lower=1> M; 				// number of areas
	int y[M];
	vector[M] E;
	matrix[M, 2] X;
	// sparse components of C matrix
	vector[M] C_eigenvalues;
	int nC_w;
	vector[nC_w] C_w; 
	int C_v[nC_w]; 
	int C_u[M+1]; 
	int offD_id_C_w[nC_w - M];		// indexes for off diagonal terms
	int D_id_C_w[M]; 				// indexes for diagonal terms - length M
}
parameters{
	// precision and sigma must be matrices
	vector<lower=0>[1] tau2_v; // marginal variance
	real<lower=0,upper=1> rho;
	vector[2] bbeta;
	vector[M] phi;
}
transformed parameters{
	// To generalize the MLCAR we must use matrices
	matrix[M,1] phi_mat = to_matrix( phi, M, 1);
	// precision is inverse of tau2_v
	matrix[1,1] Omega_R = to_matrix( 1.0 ./ tau2_v, 1, 1 ); 
	// `Sigma_R` is covariance matrix
	matrix[1,1] Sigma_R = to_matrix( tau2_v, 1, 1 );
}
model{
	// likelihood model
	y ~ poisson_log(log(E) + X * bbeta + phi);
	
	// priors
	tau2_v[1] ~ inv_gamma( 1, 0.5); // variance
	bbeta ~ std_normal();
	rho ~ uniform( 0, 1);
	
	// unit Leroux CAR prior
	target += MLCAR_lpdf( phi_mat | rho, Omega_R, Sigma_R,
	C_w, C_v, C_u,     
	offD_id_C_w, D_id_C_w, C_eigenvalues, M, 1);
}
generated quantities{
real tau2 = tau2_v[1];
}

