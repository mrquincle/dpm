% Plots the probability density function of the Pareto distribution

addpath('../inference')
addpath('../inference/prior/pareto')

alpha=5;
N=1000;
lambda_r=2;
lambda_l=-6;

% this is not what we want, it would lead to a different shape for the
% pdf on the right side versus the left side
%p1=rnd_pareto(lambda_r, alpha, 1, N);
%p2=rnd_pareto(lambda_l, alpha, 1, N);
%p=[p1 p2];

[p1 p2]=rnd_pair_two_sided_pareto(lambda_l, lambda_r, alpha, 0, N);

p=[p1 p2];

hist(p,1000);

return

x=[0:0.001:5];

k=3;
b=1;

% p(x) = k*b^k / x^(k+1) for x < b

y=(x > b) * (k*b^k) ./ x.^(k+1);

plot(x,y,'b-');

hold on

k=2;
b=2;
y=(x > b) * (k*b^k) ./ x.^(k+1);

plot(x,y,'g-');
