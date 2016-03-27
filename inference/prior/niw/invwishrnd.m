function S = invwishrnd(n,lambda)

	% Sample from an inverse Wishart distribution
	% n : degrees of freedom
	% lambda : scale parameter

	iS=wishrnd(n,lambda^-1);
	S=inv(iS);
end
