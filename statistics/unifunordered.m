% --Function: Y = unifunorderedpdf(X, A, B)
%    
%    This function is a small adaptation to the unifpdf function. That function
%    only sets Y = 1 when A <= X <= B. This function unifunorderedpdf does set
%    Y = 1 when A <= X <= B or when B <= X <= A.
%
%    Note: this probability density function is now to anymore identifiable 
%    with respect to the values of A and B. unifunorderedpdf(X, A, B) and
%    unifunorderedpdf(X, B, A) have the same probability.
%
function Y = unifunorderedpdf(X, A, B)
	C = min(A, B);
	D = max(A, B);
	Y = unifpdf(X, C, D);
end
