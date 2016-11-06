% --Function: n = unifcenteredpdf(X, MU, SIGMA)
%    
%    The function unifpdf already takes arguments A and B as vectors as long
%    as they are the same size as X. This function swaps values in A and B
%    in such way that the order doesn't matter.
%
function Y = unifcenteredpdf(X, MU, SIGMA)
	assert (SIGMA > 0);
	Y = unifpdf(X-MU, -SIGMA, +SIGMA);
end

