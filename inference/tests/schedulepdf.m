function Y = schedulepdf(X, weights0, weights1, weights2, a, b, mu1, kappa1, mu2, kappa2)

	pdf0 = unifpdf(X, a, b);
	factor0 = weights0 * pdf0;

	% X is between 0 and 24, and should be mapped to between -pi and +pi
	X = (X - 12)*pi/12;
	% same for mu1 and mu2
	mu1 = (mu1 - 12)*pi/12;
	mu2 = (mu2 - 12)*pi/12;

	factor1 = weights1 * vmpdf(X, mu1, kappa1);
	factor2 = weights2 * vmpdf(X, mu2, kappa2);

	factor = factor0 + factor1 + factor2;
	Y = factor;
end
