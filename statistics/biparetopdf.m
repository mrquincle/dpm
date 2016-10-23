% -- Function: biparetopdf(x, mu, Sigma)
%     Calculate bi-Pareto distribution
%     a > 0 is shape parameter, b = scale parameter
function result = biparetopdf(x, alpha, x_m, x_n)
	% x_m should be smaller than x_n
	if (x_m > x_n) 
		[x_m, x_n] = deal(x_n, x_m);
	end
	if (x_m == x_n)
		error('x_m should be unequal to x_n');
	end

	x_mn = (x_m + x_n) / 2;

	res0 = paretopdf(x - x_mn, alpha, x_n - x_mn);

	res1 = paretopdf(-(x - x_mn), alpha, -(x_m - x_mn));

	result = (res0 + res1)/2;

end
