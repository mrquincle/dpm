function [U E] = unifexprnd(a, b, lambda, N = 1)
	U = unifrnd(a, b, N, 1);
	%U = unifrnd(-pi, pi, N, 1);
	E = exprnd(lambda, N, 1);
end
