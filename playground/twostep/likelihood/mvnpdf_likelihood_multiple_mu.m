% -- Function File: Y = mvnpdf_likelihood_multiple_mu(X, MU, SIGMA)
%     Compute multivariate normal pdf likelihood for X given mean MU and covariance 
%     matrix SIGMA. The dimension of X is P x D, MU is P x K and SIGMA
%     is P x P x K.
%
%     See: mvnpdf_multiple_mu
%
%     Warning:
%     -- Note that the ordering in X differs from the normal mvnpdf function.

function res = mvnpdf_likelihood_multiple_mu(x, mu, sigma)
	res = sum(mvnpdf_multiple_mu(x, mu, sigma));
end
