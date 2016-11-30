% -- Function: RES = extend(S, F, N, DIR)
%     S is a struct, F is a fieldname, N is the number of copies.
%  
%     If DIR=0 we expect row-centric data, so the result will be copied in the 
%     column direction. We do not actually expect a single row, the data can 
%     also be multi-dimensional (a DxD matrix). If the data is of the form 
%     Dx(D*N) we assume the extension for N items has already been taken place. 
%     The data will not be extended in that case.
% 
%     if DIR=1 we expect column-centric data. The result will be copied in the
%     row direction. If the data is of the form (D*N)xD we assume that the
%     extension for N items has already been taken place.
%
function res = extend(S, F, N, DIR=1)
	res = extract(S, F, DIR);

	if (N < 1)
		error("extend: N should be one or larger");
	end

	[D1 D2] = size(S(1).(F));

	[R C] = size(res);

	already_right_size = false;

	switch(DIR)
	case 0
		if (R == (N * D2))
			already_right_size = true;
		end
	case 1
		if (C == (N * D2))
			already_right_size = true;
		end
	end

	switch(DIR)
	case 0
		if ~already_right_size
			res = repmat(res, N, 1);
		end
	case 1
		if ~already_right_size
			res = repmat(res, 1, N);
		end
	end
end

