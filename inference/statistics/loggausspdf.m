% -- Function: logP = loggausspdf(X, MU, SIGMA)
%     Calculate Gaussian loglikelihood for X and MU with the same dimension.
%     X is a matrix of dimensions D x N. MU has exactly the same dimensions. 
%     The individual data items X(:,i) and MU(:,i) are column-vectors and 
%     the index to pick the ith element is at the end of the list of indices.
%
%     SIGMA can be either a vector of 1 x N scalars, a D x D matrix (same for
%     each data item), a D x (D x N) matrix (different for each data item), or
%     a D x D x N nd-array (similar function as D x (D x N) matrix).
%
%     logP is a vector of 1 x N scalars
%
%     Alternate usage:
%       + MU can be a single vector of dimension D x 1 while X is D x N. In that
%       case the same MU is used for all X.
%
%     See: https://www.wikiwand.com/en/Multivariate_normal_distribution
%
%     Implementation:
%
%     Probability of the multivariate Normal distribution
%       (2*pi*Sigma)^(-1/2) * exp( (x-mu)' * inv(Sigma) * (x-mu) )^(-1/2)
%     Log-likelihood
%       -k/2 ln(2*pi) - 1/2 ln(|Sigma|) - 1/2 (x-mu)' * inv(Sigma) * (x-mu)
%
function logP = loggausspdf(x, mu, Sigma)

    % Get dimensions
    [data_D, data_N] = size(x);

    [mu_D, mu_N] = size(mu);
    
    [sigma_D, sigma_N] = size(Sigma);
    
    % Check dimension and center data
    if mu_D ~= data_D
        data_cast = 'dimension_mismatch_data_vs_parameter';
    elseif mu_N == data_N
        data_cast = 'size_match_data_vs_parameter';
    elseif mu_N == 1 
        data_cast = 'size_reshape_parameter';
    else 
        data_cast = 'size_mismatch_data_vs_parameter_number';
    end

    switch(data_cast)
    case 'dimension_mismatch_data_vs_parameter'
        error('X and MU must have the same dimension (number of rows)');
    case 'size_mismatch_data_vs_parameter_number'
        error('MU must have a single column or exactly the same number of columns as X');
    case 'size_reshape_parameter'
        x_minus_mu = bsxfun(@minus, x, mu);
    case 'size_match_data_vs_parameter'
        x_minus_mu = x - mu;
    end
    
    % Different structure to the arguments is possible
    if (sigma_D == 1)
        mode = 'unique_sigma_scalar';
    elseif (ndims(Sigma) == 2) 
        if (sigma_N == data_N * sigma_D)
            mode = 'unique_sigma_matrix_reshaped';
        elseif (sigma_N == sigma_D)
            mode = 'same_sigma_matrix';
        else 
            error('SIGMA has incorrect dimensions (%ix%i), should be either Dx(DxN) or DxD', size(Sigma));
        end
    elseif ndims(Sigma) == 3
        mode = 'unique_sigma_matrix';
    end

    switch(mode)
    case 'unique_sigma_scalar'
        if (sigma_N != data_N)
            error('Vector SIGMA should have same number of columns as X and MU');
        end
    case 'unique_sigma_matrix'
        [sigma_D, sigma_D, sigma_N] = size(Sigma);
        if (sigma_N != data_N)
            error('SIGMA should have same number of elements as X and MU');
        end
    end 

    switch(mode)
    case 'unique_sigma_scalar'
        xstand = zeros(data_N, data_D);
        R = sqrt(Sigma);
        R = repmat(R, 1, data_D);
        xstand = x_minus_mu ./ R;
        logSqrtDetSigma = log(Sigma) * data_D/2;
    case 'same_sigma_matrix'
        xstand = zeros(data_D, data_N);
        % only calculate logSqrtDetSigma once
        [xstand(:,1) logSqrtDetSigma] = precompute(x_minus_mu(:,1), Sigma);
        for i = 2:data_N
            [xstand(:,i) ~] = precompute(x_minus_mu(:,i), Sigma);
        end
    case 'unique_sigma_matrix'
        xstand = zeros(data_D, data_N);
        logSqrtDetSigma = zeros(1, data_N);
        for i = 1:data_N
            [xstand(:,i) logSqrtDetSigma(i)] = precompute(x_minus_mu(:,i), Sigma(:,:,i));
        end
    case 'unique_sigma_matrix_reshaped'
        xstand = zeros(data_D, data_N);
        logSqrtDetSigma = zeros(1, data_N);
        for i = 1:data_N
            sSigma = Sigma(:,i*2-1:i*2);
            [xstand(:,i) logSqrtDetSigma(i)] = precompute(x_minus_mu(:,i), sSigma);
        end
    end

    % Quadratic form
    xquad = sumsq(xstand, 1);

    logP = -0.5*xquad - logSqrtDetSigma - data_D * log(2*pi)/2;

    % Where Sigma is smaller than 0, set result to -Inf
    % This won't work for Sigma as matrix
    % logP(find(Sigma<=0))=-Inf;
end

% Helper function to calculate that part that is harder to formulate in one 
% stroke for an entire vector.
%
function [x_S_x logSqrtDetSigma ] = precompute(x_minus_mu, Sigma)
    [R,err] = chol(Sigma);
    if err ~= 0
        error('Sigma must be semi-definite positive');
    end
    x_S_x = x_minus_mu' * inv(R) * x_minus_mu;
    logSqrtDetSigma = sum(log(diag(R)));
end
