% Function: extend(S, F, N)
%  S is a struct, F is a fieldname
%
%  If the intermediate result is of dimension 1, we will return N copies instead.
% 
%
function res = extend(S, F, N)
	res = reshape([S.(F)],size(S(1).(F)) .* size(S))';
	if (size(res,1)==1)
		if (N ~= 1)
			res = repmat(res, N, 1);
		end
	end
end

