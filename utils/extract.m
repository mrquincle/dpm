% Function: extract(S, F)
%  S is a struct, F is a fieldname
% 
% Note: used to mix-up dimensions in case of an (:,:,k) order in which
% the field is multidimensional itself.
%
function res = extract(S, F)
	[~, K] = size(S);
	[KD DD] = size(S(1).(F));

	unordered = reshape([S.(F)]', DD, KD * K)';

	% this does not yet take into account the right order; let's correct that
	if KD > 1
		order = reshape(repmat([1:KD:K*KD]', 1, K) + [0:K-1], 1, K*KD);

		% now it's properly ordered!
		res = unordered(order, :);
	else
		res = unordered;
	end
end
