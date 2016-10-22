% --Function: n = unifpdfN(X, A, B)
%    
%    The function unifpdf already takes arguments A and B as vectors as long
%    as they are the same size as X. This function swaps values in A and B
%    in such way that the order doesn't matter.
%
function n = unifpdfN(X, A, B)
	N=length(X);

	C = min(A,B);
	D = max(A,B);
	n=unifpdf(X, C, D)';

%	for i = 1:N
%		a=A(i);
%		b=B(i);
%		x=X(i);
%		if (a > b)
%			[b,a] = deal(a,b);
%		end
%		n(i) = prod(unifpdf(x, a, b));
%	end
%	n=n';
end
