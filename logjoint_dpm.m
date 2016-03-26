function logjoint = logjoint_dpm(y, partition, alpha, hyperG0)

% LOGJOINT_DPM computes the log of the joint distribution of the partition and the data in a DPM
%  logjoint = logjoint_dpm(y, partition, alpha, hyperG0)
%--------------------------------------------------------------------------
% INPUTS
%
%   - y:        Data matrix of size p*n
%   - partition:Vector of lengths n giving the partition of the data
%   - alpha:    Scale parameter of the DP
%   - hyperG0:  Structure with the hyperparameters of the base distribution
%
% OUTPUTS
%   - logjoint: log of the joint distribution of the partition and data
%   - logprior: log of the prior on the partition
%   - logmarglik: log of the marginal likelihood
%--------------------------------------------------------------------------

n = length(partition);
[ind, ~, partition2] = unique(partition);
K = length(ind);
m = zeros(K, 1);
for j=1:length(ind)
    ind_j = find(partition==ind(j));
    m(j) = length(ind_j);
    switch(hyperG0.prior)
    case { 'NIW', 'NIG' }
        hyper(j) = hyperG0;
        for k=1:length(ind_j)
            hyper(j) = update_SS(y(:,ind_j(k)), hyper(j));
        end
    end
end

% Compute log(p(partition))
logprior = K*log(alpha) + sum(gammaln(m)) + gammaln(alpha) - gammaln(alpha + n);

% Compute log(p(y|partition))
logmarglik = 0;
switch(hyperG0.prior)
case { 'NIW', 'NIG' }
    for i=1:n
        logmarglik = logmarglik + pred(y(:,i), hyper(partition2(i)));
    end
case 'DPM_Seg'
    for i=1:n
        % TODO, we cannot just use hyperG0 here with the data
        logmarglik = logmarklik + pred(y(:,i), param(partition2(i)));
    end
end

logjoint = logprior + logmarglik;
