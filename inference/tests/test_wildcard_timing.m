function test_wildcard_timing

	N = 100;
	D = 2;
	x=rand(N, D);

	mu=rand(N, D);

	size(x)

	Sigma=[ 1 0.4; 0.5 1];


	x_minus_mu = x - mu;

	tic
	for i=1:N
		[x_S_x logSqrtDetSigma ] = precompute(x_minus_mu(:,i), Sigma);
	end
	toc


end

% Helper function to calculate that part that is harder to formulate in one 
% stroke for an entire vector.
%
function [x_S_x logSqrtDetSigma ] = precompute(x_minus_mu, Sigma)
    [R,err] = chol(Sigma);
    if err ~= 0
        error('Sigma must be semi-definite positive');
    end
    x_S_x = x_minus_mu' * inv(R) * x_minus_mu;
    logSqrtDetSigma = sum(log(diag(R)));
end

