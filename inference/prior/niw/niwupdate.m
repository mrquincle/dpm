function S = niwupdate(z, S)
	mu0 = S.mu;
	kappa0 = S.kappa;
	nu0 = S.nu;
	lambda0 = S.lambda;

	% the hyperparameters can be adapted also with multiple variables at once
	% Wishart Distributions and Inverse-Wishart Sampling [Sawyer]

	% suppose, we do not just get z, but we get a series of new samples X_1, ..., X_n, with average X
	% nu_n = nu_0 + n
	% kappa_n = kappa_0 + n
	% mu_n = (kappa_0 * mu_0) / (kappa_0 + n) + n / (kappa_0 + n) * X
	% lambda_n = lambda_0 + Q_0 + (kappa_0 * n) / (kappa_0 + n) * (X - mu_0)*(X - mu_0)'
	%
	% with Q_0 is 0 if n = 1, and more specific:
	%      Q_0 = sum_{i=1}^n (X_i - X)(X_i - X)'

	kappa1 = kappa0+1;
	nu1 = nu0+1;
	mu1 = kappa0/(kappa0+1)*mu0+1/(kappa0+1)*z;
	lambda1 = lambda0+kappa0/(kappa0+1)*(z-mu0)*(z-mu0)';
	% if isnan(lambda1)
	%     keyboard
	% end
	S.kappa = kappa1;
	S.mu = mu1;
	S.nu = nu1;
	S.lambda = lambda1;
end
