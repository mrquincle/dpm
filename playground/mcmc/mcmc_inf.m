% MCMC 
%
% MCMC can be used for (at least) two things:
% 1. Sampling from some_pdf. It forms a some_rnd generator by "spending more time" in areas in the pdf with a high
%    probability mass, then in areas with a low probability mass.
% 2. Inference, the way we are using it here.


% Problem:
%  + Data is generated and handed to us (observations). We do know the model that generated these observations, but we
%    do not know the parameters of this model (hidden variables). This model is called the likelihood function.
%  + Bayesian methodology introduces a prior on these parameters. With this additional information we can calculate
%    what the parameter values should be. A MAP estimator returns the most likely parameter values, a Bayesian
%    estimator returns probabilities for all parameter values.

% So, suppose our likelihood function is a mixture of normals. A data point is generated from a normal distribution of
% which there are several, say K.

% Now, we want to define a prior over the parameters of these normal distributions. Let the variance be the same for 
% each distribution and equal to unity. Let the mean be different for each distribution, and place them on a 
% rectangular grid as on a chess board.

% - Suppose I postulate such a structured prior. What are its hyperparameters? The distance between the means?
% - Do line segments have such a regular structure to them, or is it more irregular?
% - If I set the structure with means "exactly right" then the data should not adjust these means anymore.

% - The prior does not need to have this structure, but can be more general. 
% - However, in that case the MCMC sampler has to be smart! There are two ways in which structure can be hidden!!!
%   - Traditionally, in the prior!
%   - Computationally, in the sampler!

% Priors should be "wrong". They should not be super sophisticated... they are just a prior notion! Often they are only
% chosen for computational convenience.
% So, let's assume a 2D Cauchy prior for each mean. Each mean is considered to be independent!

clear all;
clf;

plot_cauchy = false;
plot_chisquared = false;

% Get data
load pattern100_sigma0_1.data dataset K

N = 40;

% Cauchy has values that are significantly further than the rest. It's that bad, that its mean and variance do not 
% even exist for that reason. (Such a single outlier messes up the mean in such an extreme way that the longer you
% sample from a Cauchy distribution, the mean never stabilizes: a darker swan can always appear.)
if (plot_cauchy)
	location = 0; % should be a real
	scale = 1; % should be >0
	p=cauchy_rnd(location, scale, N, 2);
	plot(p(:,1),p(:,2),'.');
%	hist(p);
end

if (plot_chisquared)
	hold on;
	dof = 4; % should be >0 and integer
	p=chi2rnd(dof, N, 2);
	plot(p(:,1),p(:,2),'r.');
end

data = cat(1, dataset.data);

% Let's assume a window, 10x10 positive quadrant, and sample uniformly from it for the means
window_min = 0;
window_max = 10;
prior.mu = unifrnd(window_min, window_max, 2, K);

% We can get the probability of the data given this mu
mu(1,:) = prior.mu;

p(1) = mvnpdfK(mu(1,:), data);
for t=1:T
	mu(t+1,:) = mu(t,:) + proposal;

	p(t+1) = mvnpdfK(testmu, data);
	
	alpha = p(t+1)/p(t);
	u = rand(0,1);
	if (alpha <= u)
		
		p(t+1) = p(t); %reject
	end
end


