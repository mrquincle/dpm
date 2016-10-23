% Function: extract(S, F)
%  S is a struct, F is a fieldname
% 
%
function res = extract(S, F)
	res = reshape([S.(F)],size(S(1).(F)) .* size(S))';
end
