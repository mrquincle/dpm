function [partition_samples, partition_estimate, similarity] = gibbsDPM(y, hyperG0, alpha, niter, type_algo, doPlot)

%
% GIBBSDPM runs a Gibbs sampler for the Dirichlet mixture of Gaussians
%   [partition_samples, partition_estimate, similarity] = gibbsDPM(y, hyperG0, alpha, niter, type_algo, doPlot)
%
%--------------------------------------------------------------------------
% INPUTS
%   - y:        Data matrix of size p*n where n is the number of data
%               points and p the dimension
%   - hyperG0:  Hyperparameter of the base distribution G0 (normal-inverse Wishart)
%               Structure with the following fields: mu, tau, nu, lambda
%   - alpha:    Scale parameter of the Dirichlet process
%   - niter:    Number of MCMC iterations. The sampler performs niter/2
%               burn-in iterations and then stores nsamples=niter/2 samples
%   - type_algo: String of characters. Type of Gibbs samplers used.
%               Possible values:
%                   - 'BMQ': Sampler based on the Blackwell-MacQueen urn
%                   - 'CRP': Sampler based on the Chinese Restaurant Process
%                   - 'collapsedCRP': Collapsed Gibbs sampler based on the Chinese Restaurant Process
%                   - 'slice': Slice sampler
%   - doplot:   0: no plot; 1: plot partition at each iteration and each
%               data; 2: plot partition at each iteration
%
% OUTPUTS
%   - partition_samples: Matrix of size n*n_samples partition_samples(k,i)
%               gives the cluster assignment of data k for sample i
%   - partition_estimate: Vector of length n giving a Bayesian point estimate of
%               the partition
%   - similarity: Matrix of size n*n. similarity(i,j) provides
%               the number of MCMC samples where i and j are in the same cluster
%--------------------------------------------------------------------------
% EXAMPLE
%
%
%--------------------------------------------------------------------------

%
% References
% For algorithms 'BMQ', 'CRP' and 'collapsedCRP', respectively algorithms
% 1, 2 and 3 in
% R.M. Neal. Markov chain sampling methods for Dirichlet process methods.
% Journal of Computational and Graphical statistics, vol.9, 2000.
%
% For the algorithm 'slice':
% M. Kalli, J. Griffin, S. Walker. Slice sampling mixture models.
% Statistics and Computing, 2011.
% M.D. Fall, E. Barat. Gibbs sampling methods for Pitman-Yor mixture
% models. Technical Report, 2012.
%--------------------------------------------------------------------------
% Copyright François Caron, University of Oxford, 2014
%--------------------------------------------------------------------------

switch(type_algo)
    case 'BMQ' % Based on Blackwell-MacQueen urn model
        partition_samples = gibbsDPM_algo1(y, hyperG0, alpha, niter, doPlot);
    case 'CRP' % Based on the Chinese Restaurant Process
        partition_samples = gibbsDPM_algo2(y, hyperG0, alpha, niter, doPlot);
    case 'collapsedCRP' % Collapsed Gibbs based on the Chinese Restaurant Process
        partition_samples = gibbsDPM_algo3(y, hyperG0, alpha, niter, doPlot);
    case 'slicesampler' % Using slice sampling
        partition_samples = gibbsDPM_algo4(y, hyperG0, alpha, niter, doPlot);
    case 'auxiliaryvars' % Auxiliary variables
        partition_samples = gibbsDPM_algo8(y, hyperG0, alpha, niter, doPlot);
    otherwise
        error('Unknown algorithm ''%s''', type_algo);
end

% Get a Bayesian point estimate based on the posterior similarity matrix
% and Binder loss function
[partition_estimate, cost, similarity] = cluster_est_binder(partition_samples);
