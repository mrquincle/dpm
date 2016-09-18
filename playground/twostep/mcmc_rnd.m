% Run MCMC to generate samples from a probability density function. 

% Here we can use the mvnpdfK function that knows for each data point x, the probability p which with it occurs.
% We now like to start with this pdf and generate data points x. In other words, MCMC is used to get the 
% function mvnrndK using mvnpdfK as input.

clear all
clear
clf

load pattern100_sigma0_1.data dataset K

limit_min=0;
limit_max=12;

data = cat(1, dataset.data);
control.mu = cat(1, dataset.mu);
control.sigma = cat(1, dataset.sigma);

printf("Size of the dataset: %i\n", length(data));

use_resets = false;
use_structure = true;

perform_inference = false;

% scale standard deviation of the proposal distribution
stddev_scale_0 = 0.01;
random_scale_0 = 10;

big_jump_probability = 0.5;

if (use_structure) 

else
	% When not using structural input, adjust stddev_scale to "better" value
	stddev_scale_0 = stddev_scale_0 * 10;
end
stddev_scale = stddev_scale_0;

T=50000;
xwalk = zeros(T+1,2);
pwalk = zeros(T+1,1);

if (perform_inference)
	muwalk = zeros(T+1,K*2);
%	sigmawalk = zeros(T+1,K*4);

	% Let's assume a window, 10x10 positive quadrant, and sample uniformly from it for the means
	window_min = 0;
	window_max = 10;
	prior.mu = unifrnd(window_min, window_max, 2, K);
end

xwalk(1,:) = rand(2,1) * random_scale_0;
pwalk(1) = mvnpdfK(xwalk(1,:), control.mu, control.sigma);

if (perform_inference)
	muwalk(1,:) = prior.mu;
end

printf("The random walk starts at: [%i %i]\n", xwalk(1,:));

modes_found = false;

Treset=T/10;
Tmodes=1000;

burnin=0.1;

% We should separate the procedure:

% 1) Get a set of a samples from different modes, the only two things that are important are:
%    a) Make sure that if the samples are generated by counting rejections, that the high count is not caused by just 
%       having nowhere to go rather than being at a highly probable location.
%    b) Make sure that samples are not generated by only one mode. Do not set the average step size of the random
%       walk too low, if that's how is searched for other modes.
%    In a 2D space it is nice to have two independent jump directions: a minimum of 3 modes, e.g. one that goes 
%    vertical and another that goes horizontal (or diagonal).  
%    It might very well be useful to reset often to get to different modes. The Markov property of an MCMC chain is
%    not really useful when searching for modes.
%
% 2) Extract the modes from the samples, by maximizing/counting, etc. 
%
% 3) Come up with a jump proposal distribution for the macro-steps. Consider e.g. also multiples of these steps or 
%    fractions.

if (use_structure)

	for t = 1:Tmodes
		% reset every so now and then
		if (mod(t, Tmodes/10) == 0) 
			xwalk(t,:) = rand(2,1) * random_scale_0;
			printf("Another random walk starts at: [%i %i]\n", xwalk(t,:));
			pwalk(t) = mvnpdfK(xwalk(t,:), control.mu, control.sigma);
		end

		% we can use a standard random generator for the proposal distribution
		proppdf = mvnrnd([0 0], eye(2)*stddev_scale);

		% propose step using proposal distribution
		if (perform_inference)
			muwalk(t+1,:) = muwalk(t,:) + proppdf;
		else
			xwalk(t+1,:) = xwalk(t,:) + proppdf;
		end

		% and we can have a "complicated" pdf we sample from
		pwalk(t+1) = mvnpdfK(xwalk(t+1,:), control.mu, control.sigma);

		% calculate when to reject a step
		if (pwalk(t) != 0)
			pfrac=pwalk(t+1)/pwalk(t);
			alpha=min(1, pfrac);
			u=unifrnd(0,1);
			if (alpha <= u) 
				xwalk(t+1,:) = xwalk(t,:);
				pwalk(t+1) = pwalk(t);
			end
		end
	end

	% get some modes
	[uniq, ~, idx] = unique(xwalk(1:Tmodes,:), "rows");
	cnts = accumarray(idx(:),1);
	umodes = [uniq cnts];
	staying_time = 8;
	modes = umodes(cnts > staying_time, :);
	lmodes = size(modes, 1);
	printf("Number of potential modes: %i\n", lmodes);
	if (lmodes == 1)
		printf("A single mode, likely all steps have been rejected and no good starting point has been found, failure!\n");
		return;
	elif (lmodes < 5)
		lsmodes = size(umodes(cnts >= 2, :), 1);
		printf("There are %i locations with two or more rejection steps. If there are many of these locations the landscape is too flat and modes are difficult to find!\n");
	end

	% sort modes on number of cnts
	max_num_modes=10;
	ff = flipud(sortrows(modes,3))
	nmodes = min(lmodes, max_num_modes);
	topff = ff(1:nmodes, 1:2);

	% get diffs between modes as jump proposal distribution
	[tm tn] = size(topff);
	topdiff = reshape(bsxfun(@minus,reshape(topff,[1 tm tn]),reshape(topff,[tm 1 tn])),[tm*tm tn]);
	modes_found = true;

	% we need now to interpolate between modes as well
	% the largest jumps are actually not the most useful ones, they are not characteristic to the regular pattern
	%topdiff = [topdiff ; topdiff / 2];

	% get rid of jumps that are too large
	smalljumps = abs(topdiff(:,1) .* topdiff(:,2)) < 5;
	topdiff = topdiff(find(smalljumps),:);

	% get rid of jumps that are too small
	smalljumps = abs(topdiff(:,1) .* topdiff(:,2)) != 0;
	topdiff = topdiff(find(smalljumps),:);

	round(topdiff)
	if(size(topdiff,1) < 5)
		printf("No success in finding modes! Run this piece of code again, perhaps let the chains walk for a bit longer.\n");
		return
	end

	% start at place with nonzero probability
	xwalk(1,:) = ff(1,1:2);
	pwalk(1) = mvnpdfK(xwalk(1,:), control.mu, control.sigma);
end

printf("Start normal MCMC\n");

for t = 1:T

	if (mod(t, T/10) == 0) 
		printf("At t=%i (out of total T=%i)\n", t, T)
	end

	% resets 
	if (use_resets)
		if (mod(t, Treset) == 0) 
			xwalk(t,:) = rand(2,1) * random_scale_0;
			pwalk(t) = mvnpdfK(xwalk(t,:), control.mu, control.sigma);
		end
	end

	if (use_structure)

		% we use ff and a standard random generator
		proppdf = mvnrnd([0 0], eye(2)*stddev_scale);
		
		% big step only once every step
		if (rand > big_jump_probability)
			% we pick one randomly in topdiff
			topr = randi(length(topdiff));

			proppdf = proppdf + topdiff(topr,:); 
		end
	else
		% we can use a standard random generator
		proppdf = mvnrnd([0 0], eye(2)*stddev_scale);
	end

	xwalk(t+1,:) = xwalk(t,:) + proppdf;

	% and we can have a "complicated" pdf we sample from
	pwalk(t+1) = mvnpdfK(xwalk(t+1,:), control.mu, control.sigma);

	if (pwalk(t) != 0)
		pfrac=pwalk(t+1)/pwalk(t);
		alpha=min(1, pfrac);
		u=unifrnd(0,1);
		if (alpha <= u) 
			% reject
			xwalk(t+1,:) = xwalk(t,:);
			pwalk(t+1) = pwalk(t);
		end
	end
end

if (use_resets)
	% burn every first items after a Treset
	selectset = find(repmat([zeros(1,burnin*Treset) ones(1,(1-burnin)*Treset)],1,T/Treset));
else
	% burn first items overall
	selectset = find([zeros(1,burnin*T) ones(1,(1-burnin)*T)]);
end
xvalues=xwalk(selectset,:,:);

% downsample
downsample=4;
xvalues=xvalues(1:downsample:end,:);

figure(1)
subplot(1,3,1)
plot(xvalues(:,1),xvalues(:,2),'.');
limits=[limit_min limit_max limit_min limit_max];
axis(limits, "square", "manual");
title("Samples");

subplot(1,3,2)

resolution=201;
t = linspace(limit_min, limit_max, resolution);
[cx cy] = meshgrid(t, t);
xx = [cx(:) cy(:)];

y = mvnpdfK(xx, control.mu, control.sigma);
yy = reshape(y, [resolution resolution]);
contour(t,t,yy);

%plot3(cx(:),cy(:)',y);
%mesh(data(:,1),data(:,2));

axis(limits, "square", "manual");
title("Probability density function");

subplot(1,3,3)
plot(data(:,1), data(:,2), '.');
axis(limits, "square", "manual");
title("Input data");

hold off;

