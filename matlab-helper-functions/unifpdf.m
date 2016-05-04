function y = unifpdf(x,a,b)
%UNIFPDF Uniform (continuous) probability density function (pdf).
%	Y = UNIFPDF(X,A,B) returns the continuous uniform pdf on the
%	interval [A,B] at the values in X. By default A = 0 and B = 1.
%	
%	The size of Y is the common size of the input arguments. A scalar input
%	functions as a constant matrix of the same size as the other inputs.	 

%	Reference:
%	   [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%	   Functions", Government Printing Office, 1964, 26.1.34.

%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:57:07 $

if nargin < 1
 error('Requires at least one input argument.'); 
end

if nargin == 1
    a = 0;
    b = 1;
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('The arguments must be the same size or be scalars.');
end

% Initialize Y to zero.
y = zeros(size(x));

k1 = find(a >= b);
if any(k1)
    y(k1) = NaN * ones(size(x));
end

k = find(x >= a & x <= b & a < b);
if any(k),
    y(k) = 1 ./ (b(k) - a(k));
end




