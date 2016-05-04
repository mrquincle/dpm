% Plots the probability density function of the Pareto distribution

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
