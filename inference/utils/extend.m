% Function: extend(S, F, N)
%  S is a struct, F is a fieldname
%
%  If the intermediate result is of dimension 1, we will return N copies instead.
% 
%
function res = extend(S, F, N, DIR=1)
	res = extract(S, F, DIR);
	R = size(res,1);

	if (R == 1)
		if (N ~= 1)
			res = repmat(res, N, 1);
		end
	else 
		if (R != N) && (size(res,2) == 1)
			warning("extend: input is expected as row-vectors");
		end
	end
end

