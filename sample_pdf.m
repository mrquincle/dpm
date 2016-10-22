% Sample probability density function
function R = sample_pdf(S)

	switch(S.prior)
	case 'NIW'
		[R.mu, R.Sigma] = normalinvwishrnd(S);
	case 'NIG'
		[R.mu, R.Sigma] = nigrnd(S);
	case 'DPM_Seg'
		[R.mu, R.Sigma] = nigrnd(S);
		[R.a, R.b] = pareto2rnd(S);
		mu_shift = normalinvwishrnd(S.shift);
		R.a = R.a + mu_shift;
		R.b = R.b + mu_shift;
	case 'Schedule'
		% Proper hyper parameters for a von-mises?
		% use a uniform distribution (could have been a von-mises with kappa = 0)
		% use an exp distr as prior for kappa (lower hyper parameter lambda is less uniform)
		[R.mu1 R.kappa1] = unifexprnd(S.a, S.b, S.lambda);
		[R.mu2 R.kappa2] = unifexprnd(S.a, S.b, S.lambda);

		% distribution that is series of 0 and 1, with each 50% change, but that never has 0 0 0 as
		% output
		% we use for that a zero-deflated Bernoulli distribution
		weights = normalize (zerodeflbernrnd(S.p, 3) .* drchrnd(S.alpha, 3, 1)' );
		R.weights0 = weights(1);
		R.weights1 = weights(2);
		R.weights2 = weights(3);
		% 
		R.a = S.a;
		R.b = S.b;

%		R.mu1 = round(R.mu1);
		
	otherwise
		error('Unknown type of prior');
	end
end
