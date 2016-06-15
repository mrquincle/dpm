% Downdate sufficient statistics S with data z
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
function R = downdate_SS(z, k, S)
    z_k = z(:,k);
	switch (S.prior)
	case 'NIW'
		R = niwdowndate(z_k, S);
	case 'NIG'
		R = nigdowndate(z_k, S);
    case 'DPM_Seg'
        T = nigdowndate(z_k, S);
        % remove observation k from sufficient statistics for Pareto
        T.ind(T.ind==k) == [];
        R = paretodowndate(z_k, T);
        %error('Sufficient statistics is assumed to be max(x_i) and when deleting z it cannot be known what the new max(x_i) should be');
	otherwise
		error('Unknown type of prior');
	end
end
