% Function to test how a variable number of arguments is handled.
% For example, run:
%
%   test_varargin([0 1], 2, 3)
% 
function test_varargin(varargin) 
	printf("Arguments:\n");
	% Returns several ans= results
	varargin{:}
	
	printf("Vector of arguments:\n");
	% This "squishes" all arguments to single vector
	vargs=[varargin{:}]

	printf("Sum all arguments:\n");
	sumall=sum([varargin{:}])

	printf("Sum all arguments except for first:\n");
	% This "squishes" all except first argument (returns 2 + 3, not 1 + 2 + 3)
	sumallexceptfirst=sum([varargin{2:end}])

	test_f("Sum all again:\n", varargin{:});

	test_f("Sum first two:\n", varargin{1:2});
end


function test_f(msg, varargin) 
	printf(msg);
	sumall=sum([varargin{:}])
end

