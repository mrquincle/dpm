% -- Function File: Y = mm_mvnpdf(X, MU, SIGMA)
%     Compute multi-modal multivariate normal pdf for each X given a vector of 
%     mean MU and a vector of covariance matrices SIGMA. The dimension of X is 
%     D x P, MU is 1 x P x K and SIGMA is P x P x K. The dimension of Y will be 
%     D x 1.
%
%     Example:
%		K = 4; P = 2; D = 10;
%		mu = rand(1,P,K);
%		sigma = repmat(eye(P),1,1,K) * 0.01;
%		x = zeros(D,P);
%		for i = 1:D
%		   k = unidrnd(K);
%		   x(i,:) = mvnrnd(mu(:,:,k), sigma(:,:,k));
%		end
%		mm_mvnpdf(x, mu, sigma);
%
%     This function is ordinarily used if there are multiple (K) multivariate 
%     normal distributions each defined through the parameters MU(:,k) and 
%     SIGMA(:,:,k). The observations X are of dimension P (just as MU), and
%     there are D of such observations. The observations are supposed to be 
%     sampled from a multi-modal probability distribution and not from 
%     individual normal distributions. The modes are all of the same strength,
%     there are no weights in this mixture model.
%
%     Of course, this function can also be called for a single observation.
%     In that case there is still an average over all K mixture components, 
%     but only for this single observation.
%
%     Implementation considerations:
%     -- Set K at the end, so sigma(:,:,i) returns continguous memory, rather
%        than sigma(i,:,:) which returns non-continguous memory.
%     -- Allocate first, assign with (i), and sum at the end at once.
%
function res = mm_mvnpdf(x, mu, sigma)
	% Check for PxK part in mu and sigma
	assert(size(mu)(2:3) == size(sigma)(2:3));
	
	% Get dimensions
	P = size(sigma, 1);
	K = size(sigma, 3);
	D = size(x, 1);

	% Uncomment to print dimensions
	% printf("Dimensionality of the data, mu, and covar: %i\n", P);
	% printf("Number of observations: %i\n", D);
	% printf("Number of modes in multimodal mvnpdf: %i\n", K);
	
	% Compare P part in X and SIGMA
	assert(size(x, 2) == P);

	% Construct result matrix
	result = zeros(D,K);

	% Calculate all K multivariate normals given all N data points
	for i = 1:K
		result(:,i) = mvnpdf(x, mu(:,:,i), sigma(:,:,i));
	end
	% Average over the K multivariate normals (average over the rows)
	res = mean(result, 2);
end
