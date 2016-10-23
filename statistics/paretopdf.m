% -- Function: paretopdf(x, mu, Sigma)
%     Calculate Pareto distribution
%     a > 0 is shape parameter, x_m = scale parameter
function result = paretopdf(x, alpha, x_m)
	result = (x > x_m) * (alpha*x_m.^alpha) ./ x.^(alpha+1);
end
