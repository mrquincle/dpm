% Test how to work with function handles
% 
% Suppose, we have a function we want to call with a function handle:
%   call_function(@dist_function, 1, 2);
% How can we without changing the signature of the call_function method return
% a sum over the results of dist_function?
%
% In the following, we create a sum_function that can exactly be used for that.
% The call becomes:
%   call_function(@sum_function, @dist_function, 1, 2);
%

function result = test_function() 
	f = @dist_function;
	y = call_function(f, [1 2], [1 2]);
	y
	
	% now how to set up the call such that call_function receives 
	% sum(dist_function([1 2], [1 2])) as argument
	h = @sum_function;
	w = call_function(h, f, [1 2], [1 2]);
	w
end

function result = sum_function(func, varargin)
	intermediate = func(varargin{:})
	result = sum(func(varargin{:}));
end

function result = dist_function(a, b) 
	result = sqrt(a.*a + b.*b);
end

function result = call_function(func, varargin)
	result = func(varargin{:});
end


