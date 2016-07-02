% Generate samples from Angular distribution
%
% -- Function: [x y] = angularrnd(S, T)
%     Return a single sample from a Gaussian distribution according to the
%     hyperparameters in S, transformed according to T.
%     * S.sigma >= 0: Gaussian noise (scalar)
%     * T.theta with 0 <= T.theta < 2*pi: angle of the line (scalar)
%     * T.p >= 0: distance of the line with the origin (scalar)
%     * T.d1: end of line segment (scalar)
%     * T.d2: end of line segment (scalar)
% -- Function: [x y] = angularrnd(S, T, N)
%     Return N samples, see angularrnd(S, T)

function [x y] = angularrnd(S, T, N)
    if ~exist('N','var')
        N=1;
    end
    % Calculate closest point to the origin.
    LP = [T.p*cos(T.theta) T.p*sin(T.theta)];

    % Generate points uniformly on (not yet translated line segment).
    di = T.d1 + (T.d2 - T.d1) * rand(N, 1);

    % Generate orthogonal deviations from Gaussian distribution
    ni = normrnd(0, S.sigma, N, 1);

    % Generate the three components that make up the points
    % XY1 = [d1*cos(theta+pi/2) d1*sin(theta+pi/2)] + LP;
    % XY2 = [d2*cos(theta+pi/2) d2*sin(theta+pi/2)] + LP; 
    % And cos(theta+pi/2) = -sin(theta), and sin(theta+pi/2) = cos(theta)
    P0 = repmat(LP, size(di));
    P1 = [di*-sin(T.theta) di*cos(T.theta)];
    P2 = [ni*cos(T.theta) ni*sin(T.theta)];

    % Sum them to get the points
    P = P0 + P1 + P2;

    % Return values
    x = P(:, 1);
    y = P(:, 2);
end
