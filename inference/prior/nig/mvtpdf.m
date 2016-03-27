% Probability density function for multivariate Student-t distribution
%
% -- Function: Y = mvtpdf(X, MU, SIGMA, NU)
%     Compute multivariate Student-t distribution for X, mean MU, and covariance
%     matrix SIGMA. The dimension of X is N X D, MU is 1 X D, SIGMA
%     is D x D and NU is a scalar. The multivariate Student-t pdf is defined as:
%
%         Y = Gamma((NU+D)/2) * 1/Gamma(NU/2) * pi^(-D/2) * |SIGMA|^(-1/2)
%           * [ 1 + 1/NU*( (X-MU)' inv(SIGMA) (X-MU) )]^(-1/2 (NU + D))

% = gammaln((nu+d)/2) - gammaln(nu/2) - (d/2)*log(pi) -1/2*logdet(nu * sigma)
%    - 1/2*(nu+d)*log(1 + [DX' inv(Sigma) DX] /nu)
%
% *References*
%
% Wikipedia https://www.wikiwand.com/en/Multivariate_t-distribution

% Note: currently just copied from log pdf of a student-t distribution.
% from pmtk3 toolkit
% https://code.google.com/p/pmtk3/source/browse/trunk/toolbox/studentT/mvtLogpdf.m?r=88

function Y = mvtpdf(X, MU, SIGMA, NU)
	% X(i,:) is i'th case
	[N d] = size(X);
	M = repmat(MU(:)', N, 1); % replicate the mean across rows
	X = X-M;
	%X
	mahal = X'*inv(SIGMA)*X;
	%mahal
	%mahal = sum((X*inv(SIGMA)).*X,2);
	logc = gammaln((NU+d)/2) - gammaln(NU/2) - 0.5*logdet(SIGMA) ...
		- (d/2)*log(NU) - (d/2)*log(pi);
	logp = logc  -(NU+d)/2*log1p(mahal/NU);

	Y = exp(logp);
end
