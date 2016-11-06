% -- Function File: [S1, S2, ..., SN] = inds2sub (DIMS, IND1, ..., INDK)
% -- Function File: [S1, S2, ..., SN] = inds2sub (DIMS, IND_VECTOR)
%     Convert an array of linear indeces to subscripts.
%
%     This function performs ind2sub consecutively for all indices IND1, ..., 
%     INDK or all indices in the vector IND_VECTOR.
%
%     Example: 
%
%         [r, c] = inds2sub([3, 3], 8, 2)
%             => r = [2, 2]
%             => c = [3, 1]
%
%     Copyright: Anne van Rossum, Crownstone B.V., 2016
%
function [varargout] = inds2sub(dims, varargin)
	% If dimensions mismatch, make error message a little more clear
	if (nargout > length(dims))
		error("There can't be more variables for indices than the number of dimensions of the data structure");
	end

	% Cast variable list and vector of indices to same data structure
	IW = [varargin{:}];

	% Construct structure to store result
	R = zeros(length(IW), nargout);

	% Iterate over all indices and call ind2sub for each
	for k=1:length(IW)
		[temp{1:nargout}] = ind2sub(dims, IW(k)) ;
		R(k,:) = [temp{:}];
	end

	% Return result as variable output arguments just like ind2sub
	for k=1:nargout
		varargout{k} = R(:, k)'; 
	end
end

