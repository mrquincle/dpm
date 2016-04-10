% This test calculates and shows the likelihoods when we move "b" one of the 
% endpoints of our segment. As long as "b" is estimated to be before an actual
% observed value of the x-coordinate the likelihood is zero. As soon as the
% observation falls within the segment it goes suddenly goes up, and then 
% slowly diminishes.
% Somehow it would be interesting to see if the fall-off when the segment is
% considered further from the point can be faster.

addpath('..');
addpath('../inference/likelihood/normal');

method='DPM_Seg';

test_mu=false;
test_b=true;

max_x_axis=20;
max_y_axis=20;

N=100;

% Default line parameters	
mu0=[5 2];

if (test_mu)
	% Different line parameters to test
	mu1=[1 2.1];
	mu2=[2 3];
	mu=[mu0' mu1' mu2'];
else
	mu=repmat(mu0, N, 1)';
end

if (test_b)
	% test b in particular, the end of the line segment
	x_inc=max_x_axis/N;
	b_test=0:x_inc:max_x_axis-x_inc;
	b_test=b_test+x_inc;
end

N=length(mu);

% create x between a0 and b0
a0=5;
b0=10;
x=rand(1,1)*a0+b0-a0;
printf("The x generated between (%i,%i) is %i\n", a0, b0, x);
% repeat this random value for all lines
x=repmat(x, 1, N);
% append with ones, so it can be multiplied with the mu's
X=[ones(1,N);x];

% assume the method calculates the optimal mu, then use this to come up with
% y, and measure if this indeed corresponds with the likelihood expectations
y=sum(X.*mu,1);

% collect total data, to be used in algorithm
data=[X; y];

clear R
for i = 1:N
	R(i).mu = mu(:,1);
	R(i).Sigma = 0.1;
	R(i).a=0;
	R(i).b=b_test(i);
end

% show the line going from [0,10] on the x-axis and having the corresponding
% y-coordinates derived via mu0

% calculate X0=0 and Xn=10
x0=0; xn=10;
X0=[1;x0]; 
Xn=[1;xn];
y0=sum(X0.*mu0',1);
yn=sum(Xn.*mu0',1);
subplot(2,1,1)
plot([x0, xn], [y0, yn]);
axis([0 20], [0 20]);
hold on;

% point x
xp=X(2,1);
yp=y(1);
plot(xp, yp, 'o');

% calculate the likelihoods with the given data
n=likelihoods(method, data, R);

% show in the second plot the likelihood with b_test
subplot(2,1,2);

plot(b_test', n, 'r');

hold off;
%b=waitforbuttonpress();
