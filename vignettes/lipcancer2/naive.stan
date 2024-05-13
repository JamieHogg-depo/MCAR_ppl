data{
	int<lower=1> M; 				// number of areas
	int y[M];
	vector[M] E;
	matrix[M, 2] X;
	matrix[M, M] C;
}
transformed data{
	vector[M] zeros = rep_vector(0, M);
}
parameters{
	real<lower=0> tau2; // marginal variance
	real<lower=0,upper=1> rho;
	vector[2] bbeta;
	vector[M] phi;
}
model{
	// likelihood model
	y ~ poisson_log(log(E) + X * bbeta + phi);
	
	// priors
	tau2 ~ inv_gamma( 1, 0.5); // marginal variance
	bbeta ~ std_normal( );
	rho ~ uniform( 0, 1);
	
	// naive implementation of leroux
	phi ~ multi_normal_prec( zeros, ( 1.0 ./ tau2 ) * add_diag( - rho * C, 1 ) );
}

