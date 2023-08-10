functions {
real ICAR1_lpdf(
    vector x,  
    vector C_w , 
    int [] C_v , 
    int [] C_u , 
    int [] offD_id_C_w ,        
    int [] D_id_C_w ,       
    vector C_eigenvalues,       
    int M                   
    ) {                 
        vector [ num_elements(C_w) ] ImrhoC;
        vector[M] A_S;
        // Multiple off-diagonal elements by rho
        ImrhoC [ offD_id_C_w ] = - C_w[ offD_id_C_w ];
        // Calculate diagonal elements of ImrhoC
        ImrhoC [ D_id_C_w ] = 1 - C_w[ D_id_C_w ];
        A_S = csr_matrix_times_vector( M, M, ImrhoC, C_v, C_u, x );
        return -0.5 * ( 
        M*log( 2 * pi() ) 
        + dot_product(x, A_S) )
		+ normal_lpdf(sum(x) | 0, 0.001 * M);
}
real ICAR2_lpdf(vector phi, int N, int[] node1, int[] node2) {
  return -0.5 * dot_self(phi[node1] - phi[node2]) 
  + normal_lpdf(sum(phi) | 0, 0.001 * N);
}
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
real ICAR3_lpdf(
    matrix y_mat,                   // input pars in matrix form M x K
    matrix Omega_R,               // Precision matrix for factors
    vector C_w , 
    int [] C_v , 
    int [] C_u , 
    int [] offD_id_C_w ,        // indexes for off diagonal terms
    int [] D_id_C_w ,         // indexes for diagonal terms - length M
    vector D_W,
    int M,                            // Number of areas
    int K) {                          // Number of factors
        vector[M] diag_term;
        matrix[K, M] A_R = Omega_R * y_mat';
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
data {
	int<lower=0> N;
	vector[N] y;
	int which;
	matrix[1,1] prec;
	
	// Components for ICAR1
	vector[N] C_eigenvalues;
	int nC_w;
	vector[nC_w] C_w; 
	int C_v[nC_w]; 
	int C_u[N+1]; 
	int offD_id_C_w[nC_w - N];		// indexes for off diagonal terms
	int D_id_C_w[N]; 				// indexes for diagonal terms - length M
	
	// Components for ICAR2
	int<lower=0> N_edges;
	int<lower=1, upper=N> node1[N_edges]; 
	int<lower=1, upper=N> node2[N_edges];
}
parameters {
	vector[N] v;
	real alpha;
	real<lower=0,upper=1> rho;
	real<lower=0> sigma; 
}
transformed parameters{
	matrix[N,1] v_mat = to_matrix(v);
}
model{
	target += std_normal_lpdf( alpha );
	target += std_normal_lpdf( sigma ); 
	// Works well
	if(which == 1){
		target += normal_lpdf( y | alpha + v , 0.5 );
		target += ICAR1_lpdf( v | C_w, C_v, C_u, offD_id_C_w, D_id_C_w, C_eigenvalues, N );
	}
	// Identical inference but a little slower
	if(which == 2){
		target += normal_lpdf( y | alpha + v , 0.5 );
		target += ICAR2_lpdf( v | N, node1, node2);
	}
	// NOT WORKING
	if(which == 3){
		target += normal_lpdf( y | alpha + v , 0.5 );
		target += ICAR3_lpdf( v_mat | prec, C_w, C_v, C_u, offD_id_C_w, D_id_C_w, C_eigenvalues, N, 1 );
		target += normal_lpdf( sum(v_mat[,1]) | 0, 0.001 * N );
	}
	// Using LCAR as ICAR equivalent for likelihood
	if(which == 4){
		target += LCAR_lpdf (y - alpha | rho, sigma, C_w, C_v, C_u, offD_id_C_w, D_id_C_w, C_eigenvalues, N );
		target += normal_lpdf( sum(y - alpha) | 0, 0.001 * N );		
	}
}






