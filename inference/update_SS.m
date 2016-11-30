% Update sufficient statistics S with new data Z
%
% -- Function: R = update_SS(Z, k, S)
%     The data can be just observations itself or in the case of a regression
%     problem Z=(X,y) with X the independent variables (can be a vector) and y
%     the dependent variable.
%     
%     Every observation should be stored as column k in matrix Z.
%
%     Currently supported S.prior:
%      * 'NIW', a Normal inverse Wishart distribution
%      * 'NIG', a Normal inverse Gamma distribution
%      * 'Par', a Pareto distribution
%      * 'DPM_Seg', a Pareto-NIG distribution
%     Returns R, the updated sufficient statistics.
%
function R = update_SS(Z, k, S)
    z = Z(:,k);
    switch(S.prior)
    case 'NIW'
	R = niwupdate(z, S);
    case 'NIG'
	R = nigupdate(z, S);
    case 'Par'
	R = paretoupate(z, S);
    case 'DPM_Seg'
        T = nigupdate(z, S);
        R = paretoupdate(z, T);
        R.ind(end+1)=k;
    otherwise
	 error('Unknown type of prior');
    end
end
