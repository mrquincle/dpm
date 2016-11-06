% -- Function File: sum_function(func, varargin)
%     Returns sum over func(varargin).
% 
%     Suppose, we have a function we want to call with a function handle:
%         call_function(@dist_function, 1, 2);
%     How can we without changing the signature of the call_function method 
%     return a sum over the results of dist_function?
%     With sum_function this becomes:
%         call_function(@sum_function, @dist_function, 1, 2);
%

function result = sum_function(func, varargin)
	result = sum(func(varargin{:}));
end
