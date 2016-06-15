% Generate samples from the wrapped Cauchy distribution.
%
% -- Function: theta = cauchyrnd(S)
%     Return a single sample from the wrapped Cauchy distribution
%     described by its sufficient statistics S. The sufficient statistics
%     consist out of:
%     * S.mu: location vector
%     * S.gamma > 0: scale factor
% -- Function: theta = cauchyrnd(S, N)
%     Return N samples from the wrapped Cauchy distribution described by
%     its sufficient statistics S.

function theta = cauchyrnd(S, N)
	if ~exist('N','var')
		N=1;
	end

    % use uniform distribution to perform inverse transform sampling
    u=unifrnd(0,1,1,N);
    % inverse of cdf (couldn't find it online, calculated from scratch) 
    theta = 2*atan(tanh(S.gamma/2)*tan(pi*(u-1/2)))+S.mu;
end
