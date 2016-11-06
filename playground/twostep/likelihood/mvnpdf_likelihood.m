% -- Function File: Y = mvnpdf_likelihood(X, MU, SIGMA)
%     Compute multivariate normal pdf likelihood for vector X given mean MU and
%     covariance matrix SIGMA. The dimension of X is D x P, MU is 1 x P and 
%     SIGMA is P x P. The dimension of Y is 1.
%
%     A mixture of normal pdf is defined as:
%
%     Y = sum_D ( SQRT 1/( (2 pi)^P |SIGMA| exp { (X-MU)' inv(SIGMA) (X-MU) } ) )
%
%     This function is ordinarily used if there are D observations with P the 
%     dimension of the data, and all generated from a single multivariate normal
%     distribution with parameters MU and SIGMA.
%
%     Implementation considerations:
%     -- Following convention of mvnpdf function

function res = mvnpdf_likelihood(x, mu, sigma)
	res = sum(mvnpdf(x, mu, sigma));
end
