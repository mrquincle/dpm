% Calculate log-likelihood from an angular pdf
%
% -- Function: ll = logangularpdf(S, T)
%     Return log-likelihood for hyperparameters in S, and transformation
%     parameters in T for the observations (x,y).
%     * S.sigma > 0: Gaussian noise (scalar)
%     * T.theta with 0 <= T.theta < 2*pi: angle of the line (scalar)
%     * T.p >= 0: distance of the line with the origin (scalar)
%     * x: coordinates
%     * y: coordinates

function ll = logangularpdf(S, T, x, y)
    N = length(y);

    % Calculate the error with respect to the line
    %err = y * sin(T.theta) + x * cos(T.theta) - T.p;
    err = y * sin(T.theta) + x * cos(T.theta) - T.p;

    % Use this as input for the Normal distribution
    %  l(i,j)=(sigma^2)^(-n/2) * exp(-1/(2*sigma^2) * err_i'*err_i);
    %  l(i,j)=sigma^-n * exp(-1/(2*sigma^2) * err_i'*err_i);
    % becomes here in log-form:
    ll=-N*log(S.sigma) + -(1/(2*S.sigma^2) * err'*err);
end

