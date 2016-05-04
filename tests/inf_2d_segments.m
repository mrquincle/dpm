% generate samples on a line and get endpoints

lb=-4;
hb=5;
sample_cnt=100;
sample_sub_cnt = [1 3 10 sample_cnt];
emp_sigma=0.5;
avg=4.5;

lim = 10;
axis_limits=[-lim, +lim, -lim, +lim];

X=unifrnd(lb, hb, 1, sample_cnt);

Y=normrnd(avg, emp_sigma, 1, sample_cnt);

P=[X; Y];

alpha=10;

endpoint_cnt=10;

fh=figure(1);

infer_normal_distribution_with_known_variance=true;

if (infer_normal_distribution_with_known_variance)

	% priors, these are now sampled from a Normal distribution
	% if they are small, they are noninformative, if they are large they will
	% dominate the data (the data will never be able to get the segment smaller)
	pA = -normrnd(0,5);
	pB = normrnd(0,5);

	plot(X, Y, '.');
	axis(axis_limits);
	% we can use a flat prior
	% f(u|y_s) \sim exp^(-1/2emp_sigma^2/n)*(u-y_s)^2

	% or we can use a normal prior density
	% f(u|y_s) \sim exp^(-1/2emp_sigma^2)*(u-y_s)^2 * exp^(-1/2s^2)*(u-m)^2
	% this is for single observation
	% we can update m and s instead of calculating it all the time
	% posterior precision: 1/s'^2 = (emp_sigma^2 + n*s^2) / (emp_sigma^2 * s^2)
	% posterior mean: m' = m/s^2 (n/emp_sigma^2) + ....

	N=sample_cnt;

	% let us first assume variance known (emp_sigma is std.dev, so emp_sigma^2)

	prior_variance = 0.4;
	prior_precision = 1/(prior_variance^2);

	prior_mean = 1
	%empirical_mean = mean(abs(Y-mean(Y)));
	empirical_mean=mean(Y);
	empirical_variance = mean((Y-mean(Y)).^2);
	empirical_precision = 1/empirical_variance;

	% https://noppa.aalto.fi/noppa/kurssi/s-114.1310/luennot/luento_12.pdf
	% http://jackman.stanford.edu/classes/BASS/ch2.pdf
	% http://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf

	% use normal conjugacy to find out that N observations with variance emp_sigma^2 and mean x'
	% is same as single observation x' with variance emp_sigma^2/N

	known_emp_sigma = emp_sigma;
	known_precision = 1/(emp_sigma^2);
	known_mean = avg;

	post_precision = (prior_precision + N*known_precision);

	post_mean = (prior_mean * prior_precision + empirical_mean * N * known_precision) / post_precision

	post_variance = 1 / post_precision;

	% note, we actually know the precision, so we should throw away this posterior precision
	% which is just an (incorrect) estimation

	post_emp_sigma = sqrt(post_variance);

end

% actually, we would like to use a prior from
% http://jakevdp.github.io/blog/2014/06/14/\
% frequentism-and-bayesianism-4-bayesian-in-python/
% p(a,beta,emp_sigma) \sim 1/emp_sigma (1 + beta^2)^(-3/2)
% no, actually, we want to have another prior, useful for rotations

% choose the priors pA and pB very small to make them noninformative
% if they are chosen large, the data will not be able to overrule the prior
pA = -.20;
pB = .20;

fg=figure(2);

for p=1:4
	subplot(2,2,p);

	ssc=sample_sub_cnt(p);

	% Get only the first dimension of P
	S=P(1:ssc);
	Z=zeros(1:ssc, 1);
	% This could have been just a min or max operation if the points were
	% not multi-dimensional. The extremum is taken first in one dimension
	% than the other.
	[A B] = extreme_pnts(S);

	A = min(A, pA);
	B = max(B, pB);

	% sample from Pareto(b, N+K) and Pareto(a, N+K)
	% with b = max(m,pb) and a = min(m,pa)
	[S0 S1] = rnd_pair_two_sided_pareto(A, B, alpha, ssc, 1, endpoint_cnt);

	plot(S, Z, 'b*');
	hold on;
	plot(S0, 0, 'r.');
	plot(S1, 0, 'g.');
	axis(axis_limits);
	hold off;
end

% but what we actually want to do is to define a prior instead of just
% sampling from a Pareto distribution

% so, let's assume we have a Gaussian for S0
% how do we update S0?

% okay, above, solves that!! :-), yeah!!


% http://math.arizona.edu/~jwatkins/f-transform.pdf
% http://research.microsoft.com/en-us/um/people/minka/papers/minka-uniform.pdf

% now, we need to incorporate the spread of points
% suppose the points are spread across a horizontal line (but not necessary 0)
% we would have a prior across the y-axis
% and we would update it using the samples
% this is independent from above, so, pretty sure, this works

% however, now if it is a line at an angle
% we can rotate all points with the angle, such that it is a horizontal line
% we identify a prior for the rotation
% however, now we have an equation for the line
% does that not conflict with the pareto prior? no... we can do the Pareto after
% mapping all points on the 1D line





