% Dirichlet Process Mixture

% The code originates from Francois Caron from the University of Oxford.
% It has been adjusted for different applications in computer vision.
% Some of the original code might still look awkward, because of this history.
%
% There are now multiple models implemented with different priors/likelihoods.
%  - the Normal inverse Wishart prior for a Normal likelihood function with
%    unknown mean and unknown covariance matrix
%  - the Normal inverse Gamma (NIG) prior for a Normal likelihood function with
%    unknown mean and variance
%  - A Dirichlet Process prior for segments.
%
% The choice which prior is used is through the variable type_prior, which
% currently has the options 'NIG' and 'NIW'.
%
% The prior is defined in hyperG0. The NIG prior is used as a prior for the
% coefficients for a line.
%   hyperG0.prior = 'NIG';
%   hyperG0.mu = [2;0];
%   hyperG0.a = 10;
%   hyperG0.b = 0.1;
%   hyperG0.Lambda = [ 1 0.5; 0.1 1];
% Likewise it is possible to define a hyperG0.prior = 'NIW' prior. This can be
% used as a prior for a Gaussian mixture model. The hyperG0.prior = 'DPM_Seg'
% can be used as a prior for line segments. A combination of Gaussian noise
% across on top of a linear relation combined with a Pareto distribution
% restricting the line to a segment only.
%
% There are different MCMC algorithms implemented. These can be tested through
% setting the variable type_algo. The options tested for all cases are
% 'BMQ' and 'CRP'. It seems 'collapsedCRP' is implemented incorrectly.
%
% To test the classification results, RandIndex and similar metrics are used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

addpath('./performance')
addpath('./inference')
addpath('./inference/prior/nig')
addpath('./inference/prior/niw')
addpath('./inference/prior/pareto')
addpath('./inference/likelihood/normal')

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if ~isOctave
	addpath('./matlab-helper-functions')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the prior to use
prior = {'NIW'; 'NIG'; 'DPM_Seg'};
% Select the prior by picking an index of 1 or higher.
type_prior=prior{3};
fprintf('Use prior ''%s''\n', type_prior);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the dataset

outfname='twolines';
data_dir='./data/lines';
data_glob=[data_dir '/*.data.txt'];

fileList = glob(data_glob);

if (length(fileList) == 0)
	printf('There are no *.data.txt files in the %s directory\n', data_dir);
	exit
end

output_dir='./output';
if ~exist(output_dir, 'dir')
	mkdir(output_dir);
end

timestamp=datestr(now,'yyyy mmm dd HH:MM:SS.FFF');
timestamp( timestamp == ' ' ) = '_';

output_inference_file=[output_dir '/' outfname '.pnts.' timestamp '.data.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some other configuration options

% Number of iterations
niter = 1000;
niter = 200;

% Type of plots (plot at every Gibbs step, or only after each point is updated)
doPlot = 2;
doPlotRaw = 0;

% Get MAP estimate
findMAP = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet Process Mixture parameters

% The most important parameter in a DPM is the scale or concentration parameter
% called alpha. With a small alpha we will have only a few clusters. With alpha
% very large we will have many clusters.
% It is also possible to perform inference over this hyperparameter and
% postulate an (improper) prior for alpha. We won't do that in this code.
alpha = 1;
%alpha = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The inference method
% The BMQ or CRP method can be used for linear regression. The BMQ is not
% converging so fast though. The collapsedCRP method or the slicesampler
% I haven't used. The auxiliaryvars method needs to be used for the DPM_Seg
% classifier.
algorithm = {'BMQ'; 'CRP'; 'collapsedCRP'; 'slicesampler'; 'auxiliaryvars'};
type_algo = algorithm{5};
type_algo = algorithm{1};
fprintf('Use algorithm ''%s''\n', type_algo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Persistent, so we can perform post-analysis on our inference
% The save command is different between matlab/octave, so replace by fprintf
% save(output_inference_file, 'output_inference_file', 'niter', 'doPlot', 'type_prior', 'alpha', 'type_algo');
%
fid=fopen(output_inference_file, 'w');
fprintf(fid, 'output_inference_file:\n%s\n', output_inference_file);
fprintf(fid, 'niter:\n%i\n', niter);
fprintf(fid, 'doPlot:\n%i\n', doPlot);
fprintf(fid, 'type_prior:\n%s\n', type_prior);
fprintf(fid, 'alpha:\n%d\n', alpha);
fprintf(fid, 'type_algo:\n%s\n', type_algo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We perform inference over a list of datasets

for f = 1:length(fileList)

	% Clear all variables when we get a new dataset, except for the ones we
	% define here: the output file name, prior type, number of iterations,
	% plotting variable, alpha, and the algorithm type
	if (isOctave)
		clear -x fileList output_inference_file fid type_prior f niter doPlot doPlotRaw alpha type_algo findMAP isOctave
	else
		clearvars -except fileList output_inference_file fid type_prior f niter doPlot doPlotRaw alpha type_algo findMAP isOctave
	end

	% Load the data from the dataset
	data_file=fileList{f,1}
	pnts=load(data_file);

	% Assume last item is label, so remove for training
	P=pnts(:,1:end-1)';

	% Prepend the X-dimension with a constant 1. This is unique for a
	% transformation to [y-intercept slope] coordinates.
	X0=ones(1,size(P,2));
	P=[X0; P];

	% Plot the raw dataset, the input for the algorithm
	if (doPlotRaw)
		if (doPlot ~= 0)
			y=P(end,:);
			X=P(2:end-1,:);

			figure('name', 'simulated data')
			plot(X, y, '.')
			xlabel('X')
			ylabel('Y')
			xlim([-1 1]*20);
			ylim([-1 1]*20);
		end
	end

	% Set the parameters of the base distribution G0
	switch(type_prior)
	case 'NIW'
		hyperG0.prior = type_prior;
		% using column-vectors
		hyperG0.mu = [0;0];
		hyperG0.kappa = 1;
		hyperG0.nu = 4;
		hyperG0.lambda = eye(2);
	case 'NIG'
		hyperG0.prior = type_prior;
		% this mean is used as prior for the coefficients for a line, 0, 0 means a horizontal line through the origin
		% [0,1] means a horizontal line through 1, [1,0] means a diagonal line through the origin (y=x)
		hyperG0.mu = [0;0];
		hyperG0.a = 10; % > 0
		hyperG0.b = 0.1;
		hyperG0.Lambda = [ 1 0.2; 0.1 1]*0.02;
		%hyperG0.Lambda = [ 1 0.5; 0.1 0.8];
	case 'DPM_Seg'
		hyperG0.prior = type_prior;
		hyperG0.mu = [2;0];
		hyperG0.a = 10; % > 0
		hyperG0.b = 0.1;
		hyperG0.Lambda = [ 1 0.5; 0.1 0.8];
%		hyperG0.Lambda = [ 1 0.2; 0.1 1]*0.02;
		% Hyper parameters for the Pareto priors
		hyperG0.p_a = -1;
		hyperG0.p_b = 1;
		hyperG0.p_alpha = 3;
		% Hyper parameter for shift on x-axis
		hyperG0.shift.mu = 0;
		hyperG0.shift.kappa = 0.05;
		hyperG0.shift.nu = 4;
		hyperG0.shift.lambda = eye(1);
	otherwise
		error('Unknown prior ''%s''', type_prior);
	end

	% Actual Gibbs sampling
	[c_st, c_est, similarity] = gibbsDPM(P, hyperG0, alpha, niter, type_algo, doPlot);

	% Plots posterior similarity matrix (data ordered with respect to point estimate)
	if (doPlot)
		[~, ind] = sort(c_est);
		figure
		imagesc(double(similarity(ind, ind)))

		% Plot point estimate of the partition
		cmap = colormap;
		figure
		for i=1:size(P, 2)
			plot(P(2,i), P(3,i), 'o', 'color', cmap(mod(10*c_est(i),63)+1,:), 'linewidth', 3);
			hold on
		end
		hold off
		drawnow()
		title('Bayesian point estimate')
	end

	% we should somehow be able to compare the assignments with the ground truth
	% c_est(i) is the index we assign to point pnts(i,end) where end is the
	% value denoting the class
	R=[c_est pnts(:,end)];
	% use standard metrics to check the performance of the assignment
	[AR,RI,MI,HI]=RandIndex(R(:,1), R(:,2))

	if (findMAP)
		% find MAP
		for i=1:size(c_st, 2)
			logjoint(i) = logjoint_dpm(P, c_st(:,i), alpha, hyperG0);
		end
		[~, indmap] = max(logjoint);
		c_map = c_st(:, indmap);
		if (doPlot)
			figure;
			plot(logjoint)
			xlabel('MCMC iteration')
			xlabel('log joint distribution')

			% Plot MAP estimate of the partition
			cmap = colormap;
			figure
			for i=1:size(P, 2)
				plot(P(2,i), P(3,i), 'o', 'color', cmap(mod(10*c_map(i),63)+1,:), 'linewidth', 3);
				hold on
			end
			title('MAP estimate')
		end

		% we should somehow be able to compare the assignments with the ground truth
		RT=[c_map pnts(:,end)];
		[AR,RI,MI,HI]=RandIndex(RT(:,1), RT(:,2))

		clen=size(unique(c_map),1);
		chist=zeros(2, clen);
		[chist(1,:) chist(2,:)] = hist(c_map, unique(c_map));

		mus=zeros(size(hyperG0.mu, 1), clen);
		for i=1:clen
			% subset of points belonging to non-empty cluster chist[i]
			Pl=P(:,c_map==chist(2,i));
			Pu=update_SS(Pl,hyperG0);
			mus(:,i)=Pu.mu';
		end
	end

	fprintf(fid, 'data_file:\n%s\n', data_file);
	fprintf(fid, 'AR:\n%f\n', AR);
	fprintf(fid, 'RI:\n%f\n', RI);
	fprintf(fid, 'MI:\n%f\n', MI);
	fprintf(fid, 'HI:\n%f\n', HI);
	if (findMAP)
		fprintf(fid, 'chist:\n%f\n', chist);
		fprintf(fid, 'mus:\n%f\n', mus);
		fprintf(fid, 'RT:\n%f\n', RT);
	end
end

fclose(fid);
