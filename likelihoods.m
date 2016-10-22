% -- Function: likelihoods(method, data, R)
%     The type of likelihood is indicated through the "method" parameter.
%     The data itself is ordered in such way that it can be compared one-on-one
%     with the entries in R. So, parameters R(1) need to be compared with
%     data(:,1).
%     Array parameters R(i).mu are reshaped to a single matrix mu(i). The same
%     for R(i).Sigma, etc.
%     If there is only a singular R(1).mu element, this element will be copied 
%     for each data point.
%
% Supported likelihood functions:
%     method='NIW' uses a Normal Inverse Wishart as prior and has a Gaussian 
%       likelihood. Given the data it will calculate the Gaussian distance 
%       between the data and (mu,Sigma) as defined in R. Note, that this method
%       seamlessly handles multiple Gaussians. R should be heterogeneous in 
%       such a case.
%     method='NIG' uses a Normal Inverse Gamma as prior and has a Gaussian
%       likelihood defined with respect to the distance between y-coordinates
%       and the predictor vector x (or design matrix X) times the searched 
%       for beta vector (line parameter vector).
%       In ridge regression Ax=b, the y coordinates are b, the x-coordinates
%       (enriched with a ones-column) are A. The variable we're after is x.
%     method='DPM_Seg' has a similar Gaussian likelihood for the line 
%       parameters, projects this line on the x-axis and defines a uniform
%       distribution between a and b.
function res = likelihoods(method, data, R, verbose=false)
	switch(method)
	case 'NIW'
		mu = reshape([R.mu],size(R(1).mu, 1), size(R(1).mu, 2) * length(R))';
		Sigma = reshape([R.Sigma],size(R(1).Sigma, 1), size(R(1).Sigma, 2) * length(R))';
		if (size(mu,1)==1)
			copies=size(data, 2);
			if (copies ~= 1)
				% expand mu, Sigma to fit data
				mu = repmat(mu, copies, 1);
				Sigma = repmat(Sigma, copies, 1);
			end
		end
		res = exp(loggausspdf(data, mu, Sigma));
	case 'NIG'
		mu = reshape([R.mu],size(R(1).mu, 1), size(R(1).mu, 2) * length(R))';
		Sigma = reshape([R.Sigma],size(R(1).Sigma, 1), size(R(1).Sigma, 2) * length(R))';
		if (size(mu,1)==1)
			copies=size(data, 2);
			if (copies ~= 1)
				% expand mu, Sigma to fit data
				mu = repmat(mu, copies, 1);
				Sigma = repmat(Sigma, copies, 1);
			end
		end
		y=data(end,:)';
		X=data(1:end-1,:);
		mu2=sum(X.*mu',1)';
		res =exp(loggausspdf(y, mu2, Sigma));
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

        % likelihood is nonzero if the point falls on the line, zero if 
        % it does not. There is no "punishment" if the point is far from
        % both end points
		n1 = unifpdfN(Xr, a, b);

		% combine the likelihood, if points do not fall on segment we
		% have n1 equal to 0, and want the result to be 0
		res =n0 .* n1;
	
	case 'Schedule'
		N = size(data, 2);
		mu1 = extend(R, 'mu1', N);
		kappa1 = extend(R, 'kappa1', N);
		mu2 = extend(R, 'mu2', N);
		kappa2 = extend(R, 'kappa2', N);
		a = extend(R, 'a', N);
		b = extend(R, 'b', N);
		weights0 = extend(R, 'weights0', N);
		weights1 = extend(R, 'weights1', N);
		weights2 = extend(R, 'weights2', N);
	
		% scale, make sure vmpdf(sched_pi(.), mu(.), kappa) is maximum
		%mu = (mu - 12)*pi/12 - pi/12;
		mu1 = (mu1 - 13)*pi/12;
		mu2 = (mu2 - 13)*pi/12;

%		K = length(R);
		sched_pi = (-pi:(2*pi)/24:pi)(1:end-1);
		prob = zeros(1,N);
		for n=1:N
			if (verbose)
				n
			end
			% we should use schedule_pdf...
			data_n = data(:,n);
			prob(n) = 1;
			J = length(data_n);
			for j=1:J
			    if (data_n(j) != 0) 
				% sum over properly weights pdfs is also a pdf (integrating to one)
				factor1 = weights0(n) * unifpdf(j-1, a(n), b(n));
				factor2 = weights1(n) * vmpdf(sched_pi(j), mu1(n), kappa1(n));
				factor3 = weights2(n) * vmpdf(sched_pi(j), mu2(n), kappa2(n));
				factor = factor1 + factor2 + factor3;
				prob(n) = prob(n) * factor.^data_n(j);
				if (verbose)
					factor1
					factor2
					factor3
					prob(n)
				end
			    end

			end

%			prob(k) = 1;
%			data(k,:)
%			J = length(data(k,:));
%			for j=1:J
%				prob(k) = prob(k) * ( weights0(k) * unifpdf(data(k,j), a(k), b(k)) + weights1(k) * vmpdf(data(k,j), mu1(k), kappa1(k)) + weights2(k) * vmpdf(data(k,j), mu2(k), kappa2(k)) );
%			end
		end
		res = prob;
	otherwise
		error('Unknown type of prior');
	end
end
