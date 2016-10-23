% -- Function: Y = zerodeflbernrnd(P,N)
%     P > 0.
%     N > 0.
%
%     A zero-deflated Bernoulli distribution. This is like a normal Bernoulli
%     distribution (a binomial distribution with one experiment), but wih zero
%     probability assigned to a vector [0 0 .. 0] as output.
%
%     The effect is strongest with N=1, at which point Y=1 with certainty. For
%     large N, the distribution becomes the normal Bernoulli distribution.
%
%     The function can be used to pick 1 or more items out of a collection.
%
function Y = zerodeflbernrnd(P, N)
	do 
		Y = binornd(1, P, N, 1);
	until (sum(Y) != 0)
end
