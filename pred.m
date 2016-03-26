% Calculate the posterior or prior predictive of (possible) observation z
% given either a prior or the data that is observed, represented by their
% sufficient statistics S.
%
% -- Function: out = pred(z, S)
%     Return a probability that z will be observed given the sufficient
%     statistics S.
%     The sufficient statistics depends on the probability density function
%     that the sufficient statistics represent. Options for S.prior:
%      * 'NIW', a Normal inverse Wishart distribution
%      * 'NIG', a Normal inverse Gamma distribution
%     This function can be used as a prior predictive by plugging in a prior.
%      - pred(z, G0)
%     It can be used as a posterior predictive by plugging in sufficient
%     statistics updated with the observed data.
%      - pred(z, updateSS(Z, ss))
%     The form is the same for conjugate distributions.
%
% See more:
%  * https://en.wikipedia.org/wiki/Posterior_predictive_distribution
function out = pred(z, S)
	switch (S.prior)
	case 'NIW'
		mu0 = S.mu;
		kappa0 = S.kappa;
		nu0 = S.nu;
		lambda0 = S.lambda;
		p = length(mu0);

		S = lambda0*(kappa0+1)/kappa0/(nu0-p+1);
		nu = nu0-p+1;
		out = (1+(z-mu0)'*S^-1*(z-mu0)/nu)^(-(nu+p)/2)*...
			exp(gammaln((nu+p)/2))/exp(gammaln((nu)/2))*...
			(det(nu*p*S))^-.5;
	case 'NIG'
		% p(z*|Z) = Student_nu(X*mu, b/a(I + X Gamma^-1 X))
		% all z are summarized through sufficient statistics, so we need only
		% p(z*|ss) and here z* consists of y* and X*
		y=z(end,:);
		X=z(1:end-1,:);
		nu = 2*S.a;
		XGX=X' * S.Lambda * X;
		out = mvtpdf(y, X'*S.mu, S.b / S.a * (eye(size(XGX)) + XGX), nu);
	otherwise
		error('Unknown type of prior');
	end

end
