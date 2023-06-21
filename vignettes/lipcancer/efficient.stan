functions{
/**
* Log probability density of the leroux conditional autoregressive (LCAR) model
* @param x vector of random effects
* @param rho spatial dependence parameter
* @param sigma standard deviation of LCAR
* @param C_w Sparse representation of C
* @param C_v Column indices for values in C
* @param C_u Row starting indices for values in C
* @param offD_id_C_w indexes for off diagonal terms
* @param D_id_C_w indexes for diagonal terms - length M
* @param C_eigenvalues eigenvalues for C
* @param M number of areas
**
@return Log probability density
*/
real LCAR_lpdf(
    vector x,               
    real rho,                   
    real sigma,              
    vector C_w , 
    int [] C_v , 
    int [] C_u , 
    int [] offD_id_C_w ,        
    int [] D_id_C_w ,       
    vector C_eigenvalues,       
    int M                   
    ) {                 
        vector[M] ldet_C;
        vector [ num_elements(C_w) ] ImrhoC;
        vector[M] A_S;
        // Multiple off-diagonal elements by rho
        ImrhoC [ offD_id_C_w ] = - rho * C_w[ offD_id_C_w ];
        // Calculate diagonal elements of ImrhoC
        ImrhoC [ D_id_C_w ] = 1 - rho * C_w[ D_id_C_w ];
        A_S = csr_matrix_times_vector( M, M, ImrhoC, C_v, C_u, x );
        ldet_C = log1m( rho * C_eigenvalues );
        return -0.5 * ( 
        M*log( 2 * pi() ) 
        - ( M * log(1/square(sigma)) + sum( ldet_C ) ) 
        + 1/square(sigma) * dot_product(x, A_S) 
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
	real<lower=0> tau2; // marginal variance
	real<lower=0,upper=1> rho;
	vector[2] bbeta;
	vector[M] phi;
}
model{
	// likelihood model
	y ~ poisson_log(log(E) + X * bbeta + phi);
	
	// priors
	tau2 ~ inv_gamma( 1, 0.5); // variance
	bbeta ~ std_normal();
	rho ~ uniform( 0, 1);
	
	// unit Leroux CAR prior
	target += LCAR_lpdf( phi | rho, sqrt(tau2),
	C_w, C_v, C_u,     
	offD_id_C_w, D_id_C_w, C_eigenvalues, M);
}

