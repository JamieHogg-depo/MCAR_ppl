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
        // vector[M] ldet_C;
        vector [ num_elements(C_w) ] ImrhoC;
        vector[M] A_S;
        // Multiple off-diagonal elements by rho
        ImrhoC [ offD_id_C_w ] = - C_w[ offD_id_C_w ];
        // Calculate diagonal elements of ImrhoC
        ImrhoC [ D_id_C_w ] = 1 - C_w[ D_id_C_w ];
        A_S = csr_matrix_times_vector( M, M, ImrhoC, C_v, C_u, x );
        // ldet_C = log1m( C_eigenvalues );
        return -0.5 * ( 
        M*log( 2 * pi() ) 
        // - ( sum( ldet_C ) ) 
        + dot_product(x, A_S) 
        );
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
        // Alternative specification for A_S
        // Omega_S as sparse matrix multiplied by columns of y_mat 
        // using crs_matrix_times_vector
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
}
data {
	int<lower=0> N;
	vector[N] y;
	int which;
	
	// Components for ICAR1
	vector[N] C_eigenvalues;
	int nC_w;
	vector[nC_w] C_w; 
	int C_v[nC_w]; 
	int C_u[N+1]; 
	int offD_id_C_w[nC_w - N];		// indexes for off diagonal terms
	int D_id_C_w[N]; 				// indexes for diagonal terms - length M
	
	// Components for ICAR2
	int<lower=0> rho_sa2_fixed;
	int<lower=0> nu_edges_sa2;
	int<lower=1, upper=N> node1_sa2[nu_edges_sa2]; 
	int<lower=1, upper=N> node2_sa2[nu_edges_sa2];
	
	// Components for ICAR1
	vector[N] C_eigenvalues;
	int nC_w;
	vector[nC_w] C_w; 
	int C_v[nC_w]; 
	int C_u[N+1]; 
	int offD_id_C_w[nC_w - N];		// indexes for off diagonal terms
	int D_id_C_w[N]; 				// indexes for diagonal terms - length M
}
transformed data{
	matrix[1,1] prec = 1.0;
}
parameters {
	vector[N] v;
	real alpha;
}
model{
	target ~ normal_lpdf( y | alpha + v , 0.5 );
	target ~ std_normal_lpdf( alpha );
	if(which == 1) target += ICAR1_lpdf( v | C_w, C_v, C_u, offD_id_C_w, D_id_C_w, C_eigenvalues, N );
	if(which == 2) target += ICAR2_lpdf( v | N, node1_sa2, node2_sa2);
	if(which == 3) target += ICAR3_lpdf( to_matrix(v) | prec, C_w, C_v, C_u, offD_id_C_w, D_id_C_w, C_eigenvalues, N, 1 );
}






