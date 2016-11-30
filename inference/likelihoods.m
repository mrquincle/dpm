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
%
%  TODO: Check if I've run multidimensional means. 

function res = likelihoods(method, data, R, verbose=false)
	switch(method)
	case 'NIW'
		N = size(data, 2);
	
		mu = extend(R, 'mu', N);
		Sigma = extend(R, 'Sigma', N);

		% A normal probability distribution on the data (X, MU)
		res = exp(loggausspdf(data, mu, Sigma));
	case 'NIG'
		N = size(data, 2);
		mu = extend(R, 'mu', N);
		Sigma = extend(R, 'Sigma', N);

		% A normal probability distribution on the data (Y, X * MU_L)
		y = data(end,:)';
		X = data(1:end-1,:);
		mu2 = sum(X.*mu',1)';
		res = exp(loggausspdf(y, mu2, Sigma));
	case 'DPM_Seg'
		N = size(data, 2);
		mu = extend(R, 'mu', N);
		Sigma = extend(R, 'Sigma', N);
		a = extend(R, 'a', N);
		b = extend(R, 'b', N);

		% A normal probability distribution on the data (Y, X * MU_L)
		y=data(end,:)';
		X=data(1:end-1,:);
		mu2=sum(X.*mu',1)';
		n0 = exp(loggausspdf(y, mu2, Sigma));

		% TODO: polar coordinates
		Xr = (y-mu(:,1)) ./ mu(:,2);

		% A uniform distribution on [a b]
		n1 = unifunorderedpdf(Xr, a, b)';

		% Result is product of Gaussian and translated Uniform distribution
		res = n0 .* n1;
	
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

		% The following makes clever use of a discrete sched_pi vector, but it
		% more readable to just call something like schedule_pdf instead
		sched_pi = (-pi:(2*pi)/24:pi)(1:end-1);
		prob = zeros(1,N);
		for n=1:N
			if (verbose)
				n
			end
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
						disp(factor1);
						disp(factor2);
						disp(factor3);
						disp(prob(n));
					end
			    end

			end
		end
		res = prob;
	otherwise
		error('Unknown type of prior');
	end
end
