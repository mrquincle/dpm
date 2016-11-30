function S0 = niwdowndate(z, Sn)
	% Copy all fields
	S0 = Sn;

	% Update sufficient statistics of Normal Wishart distribution by removing
	% data point z

	mu0 = Sn.mu;
	kappa0 = Sn.kappa;
	nu0 = Sn.nu;
	lambda0 = Sn.lambda;

	kappa1=kappa0-1;
	nu1=nu0-1;
	mu1=kappa0/(kappa0-1)*mu0-1/(kappa0-1)*z;
	lambda1=lambda0-(kappa0-1)/kappa0*(z-mu1)*(z-mu1)';
	% if isnan(lambda1)
	%     keyboard
	% end
	S0.kappa = kappa1;
	S0.mu = mu1;
	S0.nu = nu1;
	S0.lambda = lambda1;

end
