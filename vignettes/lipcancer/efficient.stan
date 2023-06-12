functions{
/**
* Log probability density of the unit-scaled Leroux CAR
* @param x_mat input pars in matrix form (M x 1)
* @param rho Spatial dependence parameter
* @param C matrix that represents spatial neighborhood 
* @param C_eigenvalues eigenvalues of C
* @param M number of areas
* @param K number of factors
**
@return Log probability density
*/
real unitLCAR_lpdf(
    matrix x_mat,               // input pars in matrix form M x K
    real rho,                   // spatial smoothing parameter
    vector C_w , 
    int [] C_v , 
    int [] C_u , 
    int [] offD_id_C_w ,        // indexes for off diagonal terms
    int [] D_id_C_w ,       	// indexes for diagonal terms - length M
    vector C_eigenvalues,       // eigenvalues for C
    int M,                      // Number of areas
    int K) {                    // Number of factors
        vector[M] ldet_C;
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
        real sq_term = sum( x_mat' .* A_S' );
        return -0.5 * ( 
        M*K*log( 2 * pi() ) 
        +  ( - K * sum( ldet_C ) ) + sq_term 
        );
}
}
data{
	int<lower=1> M; 				// number of areas
	int y[M];
	vector[M] E;
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
	real<lower=0> tau2;
	real<lower=0,upper=1> rho;
	real bbeta;
	vector[M] phi;
}
transformed parameters{
	// phi must be a column matrix
	matrix[M,1] phi_mat = to_matrix(phi, M, 1);
}
model{
	// likelihood model
	y ~ poisson_log(log(E) + bbeta + tau2 * phi);
	
	// priors
	tau2 ~ inv_gamma( 1, 0.5);
	bbeta ~ normal( 0, 10 );
	rho ~ uniform( 0, 1);
	
	// unit Leroux CAR prior
	target += unitLCAR_lpdf( phi_mat | rho, C_w, C_v, C_u, 
	offD_id_C_w, D_id_C_w, C_eigenvalues, M, 1);
}

