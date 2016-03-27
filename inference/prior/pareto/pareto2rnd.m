% Generate samples from the 2-sided symmetric Pareto distribution
%
% -- Function: [S0, S1] = pareto2rnd(S)
%     Return a pair of samples from a 2-sided symmetric Pareto distribution
%     * S.p_alpha > 0: shape parameter
%     * S.p_a > 0: scale parameter 1
%     * S.p_b > 0: scale parameter 2
%
function [S0, S1] = pareto2rnd(S)
	[S0 S1] = rnd_pair_two_sided_pareto(S.p_a, S.p_b, S.p_alpha);
end
