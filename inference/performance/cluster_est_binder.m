function [c_est, cost, similarity] = cluster_est_binder(c, thin)

%
% CLUSTER_EST_BINDER gets a point estimate of the partition using Binder loss function
%   [c_est, cost, similarity] = CLUSTER_EST_BINDER(c, thin)
%
%   INPUTS
%   - c:    L*n matrix. c(i,j) is the cluster assignment of cluster i at
%           MCMC iteration j
%   - thin: Integer. Thinning of the MCMC output (Default==1)
%
%   OUTPUTS
%   - c_est: vector of length L. Point estimate of the partition
%   - cost: vector of length N. cost(j) is the cost associated to partition
%           c(:,j)
%   - similarity: matrix of size n*n. Similarity matrix.
%
% See also COCLUSTERING
%--------------------------------------------------------------------------
% EXAMPLE:
% c_st=[1,1,2;1,2,3;2,2,1]
% [c_est, cost, similarity] = cluster_est_binder(c_st)
%

%--------------------------------------------------------------------------
% References:
% F. Caron, Y.W. Teh, T.B. Murphy. Bayesian nonparametric Plackett-Luce 
% models for the analysis of preferences for college degree programmes. 
% To appear in Annals of Applied Statistics, 2014.
% http://arxiv.org/abs/1211.5037
% D. B. Dahl. Model-Based Clustering for Expression Data via a 
% Dirichlet Process Mixture Model, in Bayesian Inference for Gene 
% Expression and Proteomics, K.-A. Do, P. Müller, M. Vannucci (Eds.),
% Cambridge University Press, 2006.
%
% Copyright Francois Caron, University of Oxford, 2014
%--------------------------------------------------------------------------


if nargin<2
    thin = 1;
end
c = c(:, 1:thin:end);

[n, N] = size(c);

similarity = zeros(n, 'int32');
cost = zeros(N, 1);
for i=1:n
    similarity(i,:) = sum(bsxfun(@eq, c(i,:), c), 2)';
    temp1 = int32(bsxfun(@eq, c(i,:), c));
    temp2 = bsxfun(@minus, N*temp1, similarity(i,:)');
    cost = cost + 1/N*sum(abs(temp2),1)';
%     pause
% Non vectorized piece of code for reference
%     for j=i+1:n
%         similarity(i,j) = 1/N*sum((c(i, :) == c(j, :)));    
%         for k=1:N
%             cost(k) = cost(k) + abs((c(i, k)==c(j,k)) - similarity(i,j) );
%         end
%     end
end

[~, ind] = min(cost);
c_est = c(:, ind);

% Label clusters by decreasing cluster size
ind1 = unique(c_est);
J = length(ind1);
clust_size = zeros(J, 1);
for i=1:J
    clust_size(i)=sum(c_est==ind1(i));
end
[~, ind] = sort(clust_size, 'descend');
c2 = c_est;
for j=1:J
    c_est(c2==ind1(ind(j))) = j;
end