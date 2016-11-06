
weights0(1) = 0.5;
weights1(1) = 0.1;
weights2(1) = 0.4;

a(1) = 0;
b(1) = 24;

mu1(1) = 6;
kappa1(1) = 1;
mu2(1) = 18;
kappa2(1) = 1;


N=1;

prob=zeros(N,1);

for n=1:N
	factor1 = weights0(n) * unifpdf(j-1, a(n), b(n));
	factor2 = weights1(n) * vmpdf(sched_pi(j), mu1(n), kappa1(n));
	factor3 = weights2(n) * vmpdf(sched_pi(j), mu2(n), kappa2(n));
	factor = factor1 + factor2 + factor3;
	prob(n) = prob(n) * factor.^data_n(j);
end

