clear all;
addpath("../utils");

K = 4;
for i = 1:K
	S(i).M = rand(3,2);
end

S.M

SM = extract(S, 'M');

SM
