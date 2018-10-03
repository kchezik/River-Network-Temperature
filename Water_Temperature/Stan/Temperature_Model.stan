data{
	int N;
	int Y;
	int year[N];
	int R;
	int site[N];
	real y[N];
	real d[N];
	real n[N];
}

parameters{
	// Local parameters
	vector<lower=0,upper=20>[Y] alpha[R];
	vector<lower=0,upper=20>[Y] A[R];
	vector<lower=0,upper=3>[Y] S[R];
	vector<lower=0.5,upper=1>[R] tau;
	vector<lower=1,upper=2>[R] tao;
	vector<lower=0,upper=5>[R] sigma;
}

transformed parameters{
	vector[N] accum;
	real annual;
	real season;
	
	for(t in 1:N){
		// Estimate annual cosine curve.
		annual = alpha[site[t],year[t]]+ A[site[t],year[t]]*cos(2*pi()*d[t]/n[t] + tau[site[t]]*pi());
		// Estimate season cosine curve.
		season = S[site[t],year[t]]*cos(2*pi()*d[t]/(n[t]/2) + tao[site[t]]*pi());
		accum[t] = normal_lpdf(y[t]| annual+ season, sigma[site[t]]);
	}
}

model{
	// Priors
	for(i in 1:R){
		alpha[i,] ~ lognormal(2,1);
		A[i,] ~ lognormal(2,1);
	}
	tau ~ normal(0.88, 0.25);
	tao ~ normal(1.5,0.25);
	sigma ~ normal(1,0.5);
	// Sum Liklihood of the data given the model parameters.
	target += sum(accum);
}
