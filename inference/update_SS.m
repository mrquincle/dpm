% Update sufficient statistics S with new data z
%
% -- Function: R = update_SS(z, k, S)
%     The data can be just observations itself or in the case of a regression
%     problem z=(X,y) with X the independent variables (can be a vector) and y
%     the dependent variable.
%     Currently supported S.prior:
%      * 'NIW', a Normal inverse Wishart distribution
%      * 'NIG', a Normal inverse Gamma distribution
%      * 'Par', a Pareto distribution
%      * 'DPM_Seg', a Pareto-NIG distribution
%     Returns R, the updated sufficient statistics.
%
function R = update_SS(z, k, S)
    z_k = z(:,k);
    switch(S.prior)
    case 'NIW'
	R = niwupdate(z_k, S);
    case 'NIG'
	R = nigupdate(z_k, S);
    case 'Par'
	R = paretoupate(z_k, S);
    case 'DPM_Seg'
        T = nigupdate(z_k, S);
        R = paretoupdate(z_k, T);
        R.ind(end+1)=k;
    otherwise
	 error('Unknown type of prior');
    end
end
