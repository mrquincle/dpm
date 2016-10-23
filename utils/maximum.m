% -- Function File: Y = maximum (X, K)
% -- Function File: [W, IW] = maximum (X, K)
% -- Function File: [W, IW] = maximum (X, K, MODE)
%     Find maximum K values in the array X.
%
%     The built-in function max does only return the maximum value. Use this
%     function to get also the 2-th till the K-th largest value.
%
%     If the optional argument MODE is given
%
function [W, IW] = maximum(X, K, MODE="range")
	[W, IW] = sort(X, 'descend');
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
