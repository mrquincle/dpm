% -- Function: Y = drchrnd(alpha, D, N)
%
%     alpha > 0
%
%     Generate N samples from a Dirichlet distribution with dispersion factor
%     alpha and dimension D. 
%     Implemented by generated Gamma distributed variables and normalizing
%     them.
%
%     For alpha approaching zero, the Dirichlet distribution picks one of the D
%     dimensions with probability approaching one. For alpha large, the 
%     Dirichlet distribution approaches a uniform distribution.
%
function Y = drchrnd(alpha, D, N)
	alpha_v = alpha(ones(1,D));

	Y = gamrnd(repmat(alpha_v, N, 1), 1, N, D);
	Y = Y ./ repmat(sum(Y,2), 1, D);
end
