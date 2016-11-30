% Function: extract(S, F)
%  S is a struct, F is a fieldname
% 
% Note: used to mix-up dimensions in case of an (:,:,k) order in which
% the field is multidimensional itself.
%
% DIR=0 assumes S.(F) to be rows. DIR=1 assumes S.(F) to be columns. The latter
% is much faster. The default is DIR=1.
%
function res = extract(S, F, DIR=1)
    [~, K] = size(S);
    [KD DD] = size(S(1).(F));

    % correctly reshaped, but not yet in order (except when KD == 1)
    if (DIR == 0)
        unordered = reshape([S.(F)]', DD, KD * K)';
    else
        unordered = [S.(F)];
    end

    % correct the order (actually not relevant in training!)
    if (KD > 1) && (DIR == 0)
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
