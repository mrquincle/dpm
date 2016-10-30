% Function: extract(S, F)
%  S is a struct, F is a fieldname
% 
% Note: used to mix-up dimensions in case of an (:,:,k) order in which
% the field is multidimensional itself.
%
function res = extract(S, F)
	[~, K] = size(S);
	[KD DD] = size(S(1).(F));

    % correctly reshaped, but not yet in order (except when KD == 1)
	unordered = reshape([S.(F)]', DD, KD * K)';

	% correct the order (actually not relevant in training!)
	if KD > 1
		% calculate how 1:K with each of length KD would be reshaped
		order = reshape(repmat([1:KD:K*KD]', 1, KD) + [0:KD-1], 1, K*KD);

		% inverse that order
		[~, inv_order] = sort(order);

		% now it's properly ordered!
		res = unordered(inv_order, :);
	else
		res = unordered;
	end
end
