% Update sufficient statistics S with new data z
%
% -- Function: R = update_SS(z, S)
%     The data can be just observations itself or in the case of a regression
%     problem z=(X,y) with X the independent variables (can be a vector) and y
%     the dependent variable.
%     Currently supported S.prior:
%      * 'NIW', a Normal inverse Wishart distribution
%      * 'NIG', a Normal inverse Gamma distribution
%     Returns R, the updated sufficient statistics.
%
function R = update_SS(z, S)
	switch(S.prior)
	case 'NIW'
		R = niwupdate(z, S);
	case 'NIG'
		R = nigupdate(z, S);
	otherwise
		error('Unknown type of prior');
	end

end
