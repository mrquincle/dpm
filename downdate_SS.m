% Downdate sufficient statistics S with new data z
%
% -- Function: R = downdate_SS(z, S)
%     The data can be just observations itself or in the case of a regression
%     problem z=(X,y) with X the independent variables (can be a vector) and y
%     the dependent variable.
%     Currently supported S.prior:
%      * 'NIW', a Normal inverse Wishart distribution
%      * 'NIG', a Normal inverse Gamma distribution
%     Returns R, the downdated sufficient statistics. If you run update_SS and
%     then downdate_SS, it should leave the argument invariant.
%
function R = downdate_SS(z, S)
	switch (S.prior)
	case 'NIW'
		R = niwdowndate(z, S);
	case 'NIG'
		R = nigdowndate(z, S);
	otherwise
		error('Unknown type of prior');
	end
end
