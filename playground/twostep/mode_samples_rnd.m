% Given a probability function, sample from those locations of high probability density.
% Currently this is implemented by restarting an MCMC a lot of times.
% The resulting walk is correlated (the steps in an MCMC usually are) and not thinned.

function res = mode_samples_rnd(Tmodes, pdf, varargin)
	Treset=Tmodes/10;
	stddev_scale=0.01;
	random_scale_0=10;
		
	% state is a to-be-generated data point
	state_init = zeros(1, 2);
	state_dim = size(state_init);
	state_walk = zeros([state_dim Tmodes+1]);
	
	% one-dimensional array of probabilities
	prob_walk = zeros(Tmodes+1);

	Tmodes
	for t = 1:Tmodes
		% reset every so now and then
		if (mod(t, Treset) == 0) 
			printf("Last state: [%i %i]\n", state_walk(:,:,t));
			state_walk(:,:,t) = rand(state_dim) * random_scale_0;
			printf("Another random walk starts at: [%i %i]\n", state_walk(:,:,t));
			prob_walk(t) = pdf(state_walk(:,:,t), varargin{:});
		end

		% we can use a standard random generator for the proposal distribution
		proppdf = mvnrnd(state_init, eye(2)*stddev_scale);

		% propose step using proposal distribution
		state_walk(:,:,t+1) = state_walk(:,:,t) + proppdf;

		% and we can have a "complicated" pdf we sample from
		prob_walk(t+1) = pdf(state_walk(:,:,t+1), varargin{:});

		% calculate when to reject a step
		if (prob_walk(t) != 0)
			pfrac=prob_walk(t+1)/prob_walk(t);
			alpha=min(1, pfrac);
			u=unifrnd(0,1);
			if (alpha <= u) 
				state_walk(:,:,t+1) = state_walk(:,:,t);
				prob_walk(t+1) = prob_walk(t);
			end
		end
	end

	res = state_walk;
end


