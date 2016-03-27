function W = wishrnd(n,lambda)

% Sample from a Wishart distribution
% n : degrees of freedom
% lambda : scale parameter

[p p2]=size(lambda);

if p~=p2
    fprintf('Error : Matrix not square\n');
end
x = chol(lambda)' * randn(p,n);
W=x*x';

