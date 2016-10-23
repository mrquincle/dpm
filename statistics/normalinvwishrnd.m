function [mu, S] = normalinvwishrnd(hyper)

	% Sample from a normal inverse Wishart distribution whose parameters are
	% given by "hyper"

	mu0 = hyper.mu;
	kappa0 = hyper.kappa;
	nu0 = hyper.nu;
	lambda0 = hyper.lambda;

	% Sample S from an inverse Wishart distribution
	S = invwishrnd(nu0,lambda0);
	% Sample mu from a normal distribution
	mu = mu0 + chol(S/kappa0)' * randn(length(mu0),1);
end
