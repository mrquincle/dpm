% Update sufficient statistics of the Normal-Inverse Wishart distribution.
%
% -- Function: Sn = niwupdate(z, S0)
%     Return updated sufficient statistics Sn when adding a single observation.
%	
% Implementation: the update to the hyperparameters when getting a series of 
% observations X_1, ..., X_n, with average X:
%   nu_n = nu_0 + n
%   kappa_n = kappa_0 + n
%   mu_n = (kappa_0 * mu_0) / (kappa_0 + n) + n / (kappa_0 + n) * X
%   Q_0 = sum_{i=1}^n (X_i - X)(X_i - X)'
%   lambda_n = lambda_0 + Q_0 + (kappa_0 * n) / (kappa_0 + n) * (X - mu_0)*(X - mu_0)'
%
% For a single observation:
%   nu_n = nu_0 + 1
%   kappa_n = kappa_0 + 1
%   mu_n = (kappa_0 * mu_0) / (kappa_0 + 1) + 1 / (kappa_0 + 1) * X
%   Q_0 = 0
%   lambda_n = lambda_0 + Q_0 + (kappa_0 * 1) / (kappa_0 + 1) * (X - mu_0)*(X - mu_0)'
%
% See: Wishart Distributions and Inverse-Wishart Sampling [Sawyer]

function Sn = niwupdate(z, S0)
	Sn = S0;
	N = size(z, 2);

	if ~(size_equal(S0.mu, z)) 
		error("niwupdate: dimension of mu and z should be equal");
	elseif (N != 1)
		error ("niwupdate: only implemented for N=1");
	end

	Sn.kappa = S0.kappa + 1;
	Sn.nu = S0.nu + 1;
	Sn.mu = (S0.mu * S0.kappa + z) / Sn.kappa;
	Sn.lambda = S0.lambda + S0.kappa / Sn.kappa * (z - S0.mu) * (z - S0.mu)';
end
