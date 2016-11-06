% Downdate according to the Normal-Inverse Gamma distribution.
%
% -- Function: S0 = nigdowndate(z, Sn)
%     Return updated sufficient statistics S0 when removing observation z from
%     sufficient statistics Sn for the Normal-Inverse Gamma distribution. The
%     result of nigdowndate(z, nigupdate(z, S0)) is S0 again.
%     The parameter z is a matrix with in each column an observation. The
%     observations consist of a pair (X,y) with X a column-vector, and y a
%     single value (the last one in the column).

function S0 = nigdowndate(z, Sn)
	% copy fields
	S0 = Sn;

	% we assume that the last row of z is y and the first rows are X
	% so each column is an observation
	y=z(end,:);
	X=z(1:end-1,:);

	%Sn.Lambda = X*X' + S0.Lambda;
	S0.Lambda = Sn.Lambda - X*X';

	%Sn.mu = inv(Sn.Lambda)*(S0.Lambda*S0.mu + X*y');
	S0.mu = inv(S0.Lambda)*(Sn.Lambda*Sn.mu - X*y');

	n = length(y');
	%Sn.a=S0.a + n/2;
	S0.a=Sn.a - n/2;

	%Sn.b=S0.b + 1/2*( y*y' + S0.mu'*S0.Lambda*S0.mu - Sn.mu'*Sn.Lambda*Sn.mu);
	S0.b=Sn.b - 1/2*( y*y' + S0.mu'*S0.Lambda*S0.mu - Sn.mu'*Sn.Lambda*Sn.mu);
end
