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

test_mu=true;
test_b=false;

max_x_axis=20;
max_y_axis=20;

% number of test cases
N=3;

% Default line parameters	
mu0 = [0 1];
b0 = 10;

if (test_mu)
    % create ground truth
    truth.mu = repmat(mu0, N, 1)';
    truth.b = repmat(b0, N, 1)';

	% Different line parameters to test
	mu1=[7.5 0.001];
	mu2=[-1500+7.5 200];
	test.mu=[mu0' mu1' mu2'];
    test.b = truth.b;
end

if (test_b)
    N = 100;
    % create ground truth
    truth.mu = repmat(mu0, N, 1)';
    truth.b = repmat(b0, N, 1)';

	% test b in particular, the end of the line segment
	x_inc=max_x_axis/N;
	test.b=0:x_inc:max_x_axis-x_inc;
	test.b=test.b+x_inc;
    test.mu = truth.mu;
end

% create x between a0 and b0
a0=5;
b0=10;
%x=rand(1,3)*a0+b0-a0;
x=[0.25 0.5 0.75]*a0+b0-a0;
%printf("The x generated between (%i,%i) is %i\n", a0, b0, x);
% repeat this random value for all lines
%x=repmat(x, 1, N);
% append with ones, so it can be multiplied with the mu's
X=[ones(1,N);x];

% assume the method calculates the optimal mu, then use this to come up with
% y, and measure if this indeed corresponds with the likelihood expectations

y=sum(X.*truth.mu,1);

% collect total data, to be used in algorithm
test.data=[X; y];

clear R
for i = 1:N
	R(i).mu = test.mu(:,i);
	R(i).Sigma = 0.1;
	%R(i).a=a0;
    %R(i).b=test.b(i);
	R(i).a=0;
    R(i).b=15;
end

% show the line going from [a0,b0] on the x-axis and having the corresponding
% y-coordinates derived for different mu
subplot(2,1,1)
for i = 1:N
    x0=a0; xn=b0;
    X0=[1;x0]; 
    Xn=[1;xn];
    y0=sum(X0.*test.mu(:,i),1);
    yn=sum(Xn.*test.mu(:,i),1);

    plot([x0, xn], [y0, yn]);
    axis([0, 20, 0, 20]);
    hold on;
end

% point x
for i = 1:N
    xp=X(2,i);
    yp=y(i);
    plot(xp, yp, 'o');
end

% calculate the likelihoods with the given data
n=likelihoods(method, test.data, R);

% show in the second plot the likelihood with b_test
subplot(2,1,2);

if (test_mu)
    n
end

if (test_b)
    axis([0, 20, 0, 20]);
    plot(b_test', n, 'r');
end

hold off;
%b=waitforbuttonpress();
