% -- Function: loggausspdf(x, mu, Sigma)
%     Calculate Gaussian loglikelihood for x and mu with the same dimension.
%     If there is only a single mu, it will be expanded to a vector the same
%     size of x. So, this compares data items with Gaussian parameters row
%     for row. It is not meant to calculate the likelihood for D=(x0,...,xN)
%     for a single Gaussian.
function logpdf = loggausspdf(x, mu, Sigma)

    [n, d] = size(x);

    [n2,d2] = size(mu);

    % Check dimension and center data
    if d2 ~= d
        error('X and mu must have the same nb of columns');
    elseif n2 == n
        x0 = x - mu;
    elseif n2 == 1 % mean is a single row, rep it out to match data
        x0 = bsxfun(@minus,x,mu);
    else % sizes don't match
        error('mu must has 1 row or the same nb as x');
    end

    [n3,d3] = size(Sigma);
    if (d3==1)
%        if (n3==1)
            % single sigma value
 %           R = sqrt(2)*Sigma
  %          xstand = x0 / R;
   %         d=1;
    %        logSqrtDetSigma=log(Sigma);
        if (n3==n)
            % row of sigma values, different for each x and mu combination
            % for n=2
            % mvnpdf(x0,mu0,diag([sigma0,sigma0]))
            % exp((x0-mu0)*(x0-mu0)'/(-2*sigma0))*1/(sqrt((2*pi*sigma0)^2))
            % logL = ((x0-mu0)*(x0-mu0)'/(-2*sigma0)) - log(sqrt((2*pi*sigma0)^d))
            %   ((x0-mu0)*(x0-mu0)'/(-2*sigma0)) - log(2*pi*sigma0)*d/2
            %   ((x0-mu0)*(x0-mu0)'/(-2*sigma0)) - log(2*pi)*d/2 - log(sigma0)*d/2
            %   -0.5*((x0-mu0)*(x0-mu0)'/sigma0) - log(2*pi)*d/2 - log(sigma0)*d/2
            xstand = zeros(n,d);
            R = sqrt(Sigma);
            R=repmat(R,1,d);
            xstand = x0 ./ R;
            logSqrtDetSigma=log(Sigma)*d/2;
        else
            error('Vector of sigma should have same number of rows as x and mu');
        end
    elseif ndims(Sigma)==2
        [R,err] = chol(Sigma);
        % Standardized data
        xstand =  x0 / R;
        if err~= 0
            error('Sigma must be semi-definite positive');
        end
        % log(sqrt(Sigma))
        logSqrtDetSigma = sum(log(diag(R)));

    elseif ndims(Sigma)==3
        xstand = zeros(n,d);
        logSqrtDetSigma = zeros(n,1);
        for i = 1:n
            [R,err] = chol(Sigma(:,:,i));
            if err ~= 0
                error('Sigma must be semi-definite positive');
            end
            xstand(i,:) = x0(i,:) / R;
            logSqrtDetSigma(i) = sum(log(diag(R)));
        end
    end

    % Quadratic form
    xquad = sum(xstand.^2, 2);
    logpdf = -0.5*xquad - logSqrtDetSigma - d*log(2*pi)/2;

    % Where Sigma is smaller than 0, set result to -Inf
    logpdf(find(Sigma<=0))=-Inf;
end
