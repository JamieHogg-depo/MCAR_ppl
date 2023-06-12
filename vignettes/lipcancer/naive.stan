data{
	int<lower=1> M; 				// number of areas
	int y[M];
	vector[M] E;
	matrix[M, M] C;
}
transformed data{
	vector[M] zeros = rep_vector(0, M);
}
parameters{
	real<lower=0> tau2;
	real<lower=0,upper=1> rho;
	real bbeta;
	vector[M] phi;
}
model{
	// likelihood model
	y ~ poisson_log(log(E) + bbeta + phi);
	
	// priors
	tau2 ~ inv_gamma( 1, 0.5);
	bbeta ~ normal( 0, 10 );
	rho ~ uniform( 0, 1);
	
	// unit Leroux CAR prior
	phi ~ multi_normal_prec( zeros, tau2 * add_diag( - rho * C, 1 ) );
}

