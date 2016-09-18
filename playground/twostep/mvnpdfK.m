function res = mvnpdfK(x, mu, sigma)
	res = mvnpdf(x, mu(1,:), sigma(1:2,:));
	[K D] = size(mu);
	Kcheck = size(sigma, 1);
	% the 
	assert(Kcheck == D*K);
	for i = 2:K
		res = res + mvnpdf(x, mu(i,:), sigma(i*2-1:i*2,:));
	end
end
