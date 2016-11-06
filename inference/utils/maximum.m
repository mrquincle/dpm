% -- Function File: Y = maximum (X, K)
% -- Function File: [W, IW] = maximum (X, K)
% -- Function File: [W, IW] = maximum (X, K, MODE)
%     Find maximum K values in the array X.
%
%     The built-in function max does only return the maximum value. Use this
%     function to get also the 2-th till the K-th largest value.
%
%     If the optional argument MODE is given it can either take the default
%     value "range" in which maximum values 1:K are returned. Or the value
%     "singular" in which only the K-th largest value is returned.
%
%     If the data structure X is a matrix, IW will be a linear index to the
%     maximum, such that X(IW) = W.
%
%     Example: 
%       > [w, iw] = maximum([4 0 10], 2)
%       w  = [10; 4]
%       iw = [3; 1]
%
%       > [w, iw] = maximum([4 0 10], 2, "singular")
%       w  = 4
%       iw = 1
%
%       > [w, iw] = maximum([4 0; 10 1], 2)
%       w  = [10; 4]
%       iw = [2; 1] 
%
function [W, IW] = maximum(X, K, MODE="range")
	[W, IW] = sort(X(:), 'descend');
	switch(MODE)
	case "range"
		W = W(1:K);
		IW = IW(1:K);
	case "singular"
		W = W(K);
		IW = IW(K);
	otherwise
		error("Unknown MODE: %s\n", MODE);
	end
end
