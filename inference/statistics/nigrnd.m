% Generate samples from the Normal-Inverse Gamma distribution.
%
% -- Function: [mu, sigma] = nigrnd(S)
%     Return a single sample from the Normal-Inverse Gamma distribution
%     described by its sufficient statistics S. The sufficient statistics
%     consist out of:
%     * S.a > 0: shape parameter
%     * S.b > 0: scale parameter
%     * S.mu: location vector
%     * S.Lambda: precision matrix, inverse effect on spread, inverse of the
%       covariance matrix
% -- Function: [mu, sigma] = nigrnd(S, N)
%     Return N samples from the Normal-Inverse Gamma distribution described by
%     its sufficient statistics S.
%
% For more information, see:
%  * https://www.wikiwand.com/en/Normal-inverse-gamma_distribution
%  * https://en.wikipedia.org/wiki/Bayesian_linear_regression
%
% The implementation for the generation of these random variates:
%  * Sample \sigma^2 from an inverse gamma distribution with parameters a and b
%  * Sample x from a normal distribution with mean \mu and variance \sigma^2/\lambda

function [mu, sigma] = nigrnd(S, N)
	if ~exist('N','var')
		N=1;
	end
	% sigma^2 \sim InverseGamma(a,b)
	sigma2=1./gamrnd(S.a, S.b, N, 1);
	% p(beta, sigma^2) \sim N(mu, sigma^2*inv(Lambda))
	mu=zeros(size(sigma2,1),2);
	for i=1:N
		mu(i,:)=mvnrnd(S.mu, sigma2(i)*inv(S.Lambda));
	end
	mu=mu';
	sigma=sqrt(sigma2);
end
