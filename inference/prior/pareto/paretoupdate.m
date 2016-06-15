% Update sufficient statistics of the Uniform-Pareto distribution.
%
% -- Function: Sn = paretoupdate(z, S0)
%     Return updated sufficient statistics Sn when adding observation(s) z from
%     sufficient statistics S0 for the Uniform-Pareto distribution. It is 
%     impossible to downdate.
%     The parameter z is a matrix with in each column an observation. The
%     observations consist of a pair (X,y) with X a column-vector, and y a
%     single value (the last one in the column).
%

% Does not work like this, z should be two endpoints, how else to update
% both sides of the distribution
function Sn = paretoupdate(z, S0)
    % copy (all) fields
    Sn = S0;

    % we assume that the last row of z is y and the first rows are X
    % so each column is an observation
    y=z(end,:);
    X=z(1:end-1,:);

    n = length(y');

    Sn.par.alpha = S0.par.alpha + n;
    % assumes X as two coordinates, second coordinate is x-axis
    Sn.beta = max(S0.beta, X(2)); 
end
