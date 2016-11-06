% Test the mm_mvnpdf function

clf
clear all
clear

addpath("pdf");

K = 4; P = 2; D = 1000;

% Create K Gaussians at [1,1], [2,2], [3,3], ... [K,K]
%mu = repmat(1:K,2,1);

mu = rand(1,P,K) * 10;
sigma = repmat(eye(P),1,1,K) * 0.01;
x = zeros(D,P);

% Overwrite with random values
x = rand(D,P) * 10;

% Overwrite first half with Gaussians
for i = 1:D/2
   k = unidrnd(K);
   % Randomly sample from one of the Gaussians
   x(i,:) = mvnrnd(mu(:,:,k), sigma(:,:,k));
end

tic
% Calculate the probability density for each sample
result = mm_mvnpdf(x, mu, sigma);
toc

% Expect: K circles at locations of Gaussians if threshold is not really low
% Only samples that have a low (unnormalized) probability will be visualized
j=(result < 0.1);
plot(x(j,1), x(j,2), 'ob')

hold on
% Expect: K disks at locations of Gaussians
% Only samples that have a high (unnormalized) probability will be visualized
j=(result >= 1);

plot(x(j,1), x(j,2), '.r')

mu
