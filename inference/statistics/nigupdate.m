% Update sufficient statistics of the Normal-Inverse Gamma distribution.
%
% -- Function: Sn = nigupdate(z, S0)
%     Return updated sufficient statistics Sn when adding observation(s) z from
%     sufficient statistics S0 for the Normal-Inverse Gamma distribution. The
%     result of nigdowndate(z, nigupdate(z, S0)) is S0 again.
%     The parameter z is a matrix with in each column an observation. The
%     observations consist of a pair (X,y) with X a column-vector, and y a
%     single value (the last one in the column).
%
% For more information, see:
%  * Bayesian Linear Model: Gory Details
%  * Pubh7440 Notes By Sudipto Banerjee
%
% If we consider V^{-1} as a mistake: Vb^{-1} and use Lambda_n = Vb^{-1} the
% implementation is as follows:
%    p(beta, sigma^2) \propto (1/sigma^2)^{a+(n+p)/2+1}
%      \exp { -1/sigma^2 [ b* + 1/2 (beta - mu*)^T V*^{-1} (beta - mu*) ] }

function Sn = nigupdate(z, S0)
	% copy fields
	Sn = S0;

	% we assume that the last row of z is y and the first rows are X
	% so each column is an observation
	y=z(end,:);
	X=z(1:end-1,:);

	Sn.Lambda = (X*X' + S0.Lambda);
	Sn.mu = inv(Sn.Lambda)*(S0.Lambda*S0.mu + X*y');
	n = length(y');
	Sn.a=S0.a + n/2;
	Sn.b=S0.b + 1/2*( y*y' + S0.mu'*S0.Lambda*S0.mu - Sn.mu'*Sn.Lambda*Sn.mu);
end
