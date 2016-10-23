% Downdate sufficient statistics of the Uniform-Pareto distribution.
%
% -- Function: Sn = paretoupdate(z, r, S0)
%     Return sufficient statistics Sn when removing observation(s) z from
%     sufficient statistics S0 for the Uniform-Pareto distribution. 
%     The parameter z is a matrix with in each column an observation. The
%     observations consist of a pair (X,y) with X a column-vector, and y a
%     single value (the last one in the column).
%

function Sn = paretodowndate(z, r, S0)
    % copy (all) fields
    Sn = S0;

    % we assume that the last row of z is y and the first rows are X
    % so each column is an observation
    z_y=z(end,:);
    n = length(z_y');

    % remaining data items are used to calculate max again
    r_X=r(1:end-1,:);


    Sn.alpha = S0.alpha - n;
    % assumes X as two coordinates, second coordinate is x-axis
    Sn.beta = max(S0.beta, r_X(2)); 
end
