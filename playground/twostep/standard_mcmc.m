% Run MCMC on the data

clear all
clear
clf

load normal.data dataset K

data = dataset(1).data;
mu = dataset(1).mu;
sigma = dataset(1).sigma;
for i = 2:K
	data = [data; dataset(i).data];
	mu = [mu; dataset(i).mu];
	sigma = [sigma; dataset(i).sigma];
end

%plot(data(:,1),data(:,2),'o');

%y = zeros(101*101,1);
%for i = 1:K
%	sel = [i*2-1:i*2];
%	y = y + mvnpdf(x, dataset(i).mu(1,:), dataset(i).sigma(1:2,:));
%end
%	y = mvnpdfK(x, dataset, K);

% Scale has to be picked right, if too low, it takes really long to mix!
scale = 0.1;
scale = 0.01;

T=10000;
xwalk = zeros(T+1,2);
pwalk = zeros(T+1,2);

xwalk(1,:) = rand(2,1);
pwalk(1,:) = mvnpdfK(xwalk(1,:), dataset, K);

for t = 1:T
	% we can use a standard random generator
	proppdf = mvnrnd([0 0], eye(2)*scale);

	xwalk(t+1,:) = xwalk(t,:) + proppdf;

	% and we can have a "complicated" pdf we sample from
	pwalk(t+1,:) = mvnpdfK(xwalk(t+1,:), dataset, K);

	alpha=min(1, pwalk(t+1,:)/pwalk(t,:));

	u=unifrnd(0,1);
	if (alpha > u) 
		% accept
	else
		% reject
		xwalk(t+1,:) = xwalk(t,:);
		pwalk(t+1,:) = pwalk(t,:);
	end
end

% burn
burnin=T/10;
xwalk=xwalk(burnin:end,:,:);
% downsample
downsample=4;
xwalk=xwalk(1:downsample:end,:);

figure(1)
subplot(1,2,1)
plot(xwalk(:,1),xwalk(:,2),'o');
limits=[-1 1 -1 1];
axis(limits, "square", "manual");

%figure(2)
subplot(1,2,2)

t = linspace(-1, 1, 101);
[cx cy] = meshgrid(t, t);
x = [cx(:) cy(:)];

y = mvnpdfK(x, dataset, K);
yy = reshape(y, [101 101]);
contour(t,t,yy);

axis(limits, "square", "manual");

%plot3(cx(:),cy(:)',y);
%mesh(data(:,1),data(:,2));


