% Generate several normal distributions

clear all;
clear;
clf;

N=20;
sqrt_K=10;
K=sqrt_K^2;

plot_data = true;
generate_pattern = true;

if (generate_pattern)
	mu_0 = rand(1,2);
end

for i=1:K
%	N = ceil(rand*100);

	if (generate_pattern)
		[Ki, Kj] = ind2sub([sqrt_K sqrt_K], 1:K);
		mu_r = mu_0 + [Ki(i) Kj(i)];
		% the smaller sigma, the more difficult!
		sigma_r1dim = [1 1] / 10;
	else
		mu_r = rand(1,2);
		sigma_r1dim = rand(2,1);
	end

%	mu = repmat(mu_r, N, 1);

	scale=0.001;
	sigma_r = diag(sigma_r1dim * scale);
	
	data = mvnrnd(mu_r,sigma_r,N);

	if (plot_data)
		plot(data(:,1),data(:,2),'o');
		hold on;
	end
	dataset(i).N = N;
	dataset(i).mu = mu_r;
	dataset(i).sigma = sigma_r;
	dataset(i).data = data;
end

if (plot_data)
	hold off;
end

if (generate_pattern)
	save pattern100_sigma0_1.data dataset K
else
	save normal.data dataset K
end

