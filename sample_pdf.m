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
	otherwise
		error('Unknown type of prior');
	end
end
