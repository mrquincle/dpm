% Pair sampling from a two-sided Pareto distribution
%
% -- Function: [S0 S1] = rnd_two_sided_pareto(x_l, x_m, alpha)
%     Return a sample pair from a two_sided Pareto distribution using inverse
%     transform sampling.
% -- Function: [S0 S1] = rnd_two_sided_pareto(x_l, x_m, alpha, N)
%     Return a sample pair from a two_sided Pareto distribution using inverse
%     transform sampling in the context of a joint distribution for {x_l, x_m}
%     and a dataset of N observations.
% -- Function: [S0 S1] = rnd_two_sided_pareto(x_l, x_m, alpha, N, R, C)
%     Return sample pairs from a two-sided Pareto distribution in R rows and C
%     columns.
%
%     Scale parameters: x_l and x_m
%     Shape parameter: alpha > 0.
%     Number of observations: N >= 0
%
%     The two-sided Pareto distribution can be used as a prior for a uniform
%     distribution over a segment.
%
%     The symmetric Pareto distribution can be sampled from by choosing
%     x_l = -x_m. For other values the distribution will be symmetric as well,
%     but not centered around the origin.
%
%     There is only one alpha parameter. If implemented just with x_l and x_m
%     separately, a single alpha leads to a different density around x_l
%     compared to x_m. Hence, the distribution is first shifted to be around
%     the origin. Then a symmetric Pareto distribution with |x_m - x_l|/2 is
%     sampled from and the result is shifted back to the [x_l, x_m] interval.
%
%     The ordering of x_l and x_m does not matter. Compared to the normal
%     Pareto distribution they neither have to be larger than zero.
%
%     Alpha should be larger than 0, or more precisely, alpha+N should be
%     larger than 0. The prior alpha functions as a pseudo-count.
%
%     If a single point is used as input, this is considered as the right
%     extreme for S0, and the left extreme for S1. The support of this
%     distribution is the entire real line.

function [S0 S1] = rnd_pair_two_sided_pareto(x_l, x_m, alpha, N, R, C)
	if ~exist('R','var')
		R=1;
	end
	if ~exist('C','var')
		C=1;
	end
	if ~exist('N','var')
		N=0;
	end
	if (x_l == x_m)
		x = abs(x_m);
		shift0 = x_m + x; % only shift if x_m is positive
		shift1 = x_m - x; % only shift if x_m is negative
	else
		x = (x_m - x_l) / 2;
		shift0 = (x_m + x_l) / 2;
		shift1 = shift0;
	end

	S0 = rnd_pareto(-x, alpha + N, R, C) + shift0;
	S1 = rnd_pareto(x, alpha + N, R, C) + shift1;
end

