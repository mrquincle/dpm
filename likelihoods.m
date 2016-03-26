% -- Function: likelihoods(method, data, R)
%     The type of likelihood is indicated through method
%     The data itself is ordered in such way that it can be compared one-on-one
%     with the entries in R. So, parameters R(1) need to be compared with
%     data(:,1)
function n = likelihoods(method, data, R)
	switch(method)
	case 'NIW'
		mu = reshape([R.mu],size(R(1).mu, 1), size(R(1).mu, 2) * length(R))';
		Sigma = reshape([R.Sigma],size(R(1).Sigma, 1), size(R(1).Sigma, 2) * length(R))';
		if (size(mu,1)==1)
			copies=size(data, 2);
			if (copies ~= 1)
				% expand mu, Sigma, a, and b to fit data
				mu = repmat(mu, copies, 1);
				Sigma = repmat(Sigma, copies, 1);
			end
		end
		n = exp(loggausspdf(data, mu, Sigma));
	case 'NIG'
		mu = reshape([R.mu],size(R(1).mu, 1), size(R(1).mu, 2) * length(R))';
		Sigma = reshape([R.Sigma],size(R(1).Sigma, 1), size(R(1).Sigma, 2) * length(R))';
		if (size(mu,1)==1)
			copies=size(data, 2);
			if (copies ~= 1)
				% expand mu, Sigma, a, and b to fit data
				mu = repmat(mu, copies, 1);
				Sigma = repmat(Sigma, copies, 1);
			end
		end
		y=data(end,:)';
		X=data(1:end-1,:);
		mu2=sum(X.*mu',1)';
		n = exp(loggausspdf(y, mu2, Sigma));
	case 'DPM_Seg'
		mu = reshape([R.mu],size(R(1).mu, 1), size(R(1).mu, 2) * length(R))';
		Sigma = reshape([R.Sigma],size(R(1).Sigma, 1), size(R(1).Sigma, 2) * length(R))';
		a = reshape([R.a],size(R(1).a, 1), size(R(1).a, 2) * length(R))';
		b = reshape([R.b],size(R(1).b, 1), size(R(1).b, 2) * length(R))';

		% if mu is size==1
		if (size(mu,1)==1)
			copies=size(data, 2);
			if (copies ~= 1)
				% expand mu, Sigma, a, and b to fit data
				mu = repmat(mu, copies, 1);
				Sigma = repmat(Sigma, copies, 1);
				a = repmat(a, copies, 1);
				b = repmat(b, copies, 1);
			end
		end

		y=data(end,:)';
		X=data(1:end-1,:);
		mu2=sum(X.*mu',1)';
		n0 = exp(loggausspdf(y, mu2, Sigma));

		% we have to expand mu and Sigma to a and b as well

		% the line parameters correspond with mu, so we can map y to X
		% through mu as well, then we only need the second coordinate
		% to get x
		%
		% let's keep it two-dimensional for now, y=mu(1)+x*mu(2)
		% hence x=(y-mu(1)) / mu(2), thus Xr are now the data points
		% projected on the x-axis
		%
		% Note, that when the line segment is almost vertical, we get
		% in trouble
		% TODO: use polar coordinates, rotate the line rather than
		% just projecting on the x-axis
		Xr = (y-mu(:,1)) ./ mu(:,2);

		% so we have a point cloud on the horizontal axis and can use
		% the uniform likelihood to establish the contribution from the
		% points being on the segment
		% we swap a and b if a is larger than b, because we don't care
		% which endpoint is which
		n1 = unifpdfN(Xr, a, b);
		% combine the likelihood, if points do not fall on segment we
		% have n1 equal to 0, and want the result to be 0
		n = n0 .* n1;
	otherwise
		error('Unknown type of prior');
	end
end
