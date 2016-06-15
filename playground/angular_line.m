
% line parameters in normal form

% angle w.r.t. origin
theta=rand*2*pi;
theta=pi/2;
% distance with origin
p=0;
% noise
sigma=0.02;

fig1=false;

% one point is defined in particular, the closest point to the origin
LP = [p*cos(theta) p*sin(theta)];
%LP
d1=-5; % travel d1 over line from LP, this is endpoint 1
d2=-9; % travel d2 over line from LP, this is endpoint 2

% should make cos(theta+pi/2) = -sin(theta), and sin(theta+pi/2) = cos(theta), but can't be bothered
%LP1= [d1*cos(theta+pi/2) d1*sin(theta+pi/2)] + LP;
%LP2= [d2*cos(theta+pi/2) d2*sin(theta+pi/2)] + LP;

if (d2 < d1) 
	[d1, d2] = deal(d2, d1);
end
di = unifrnd(d1,d2, 100,1);
XY= [di*cos(theta+pi/2) di*sin(theta+pi/2)] + repmat(LP,size(di));
x = XY(:,1);
y = XY(:,2);

theta
p
sigma

% we can distribute/generate points according to alpha and p
% x*cos(theta) + y*sin(theta) - p = 0

%x = rand(1,100)*2-1;
%y = p/sin(theta) - x/tan(theta);
%y = (p - x*cos(theta)) / sin(theta);

if (fig1)
	plot(x, y, '.')

	hold on
	axis("equal");

	% plot LP
	plot(LP(1), LP(2), '*');
	plot(LP1(1), LP1(2), 'o');
	plot(LP2(1), LP2(2), 'o');
end

x_t = -1:0.1:1;
theta_t = theta+pi/2;
theta_t
y_t = - x_t/tan(theta_t);

%plot(x_t,y_t);

t = linspace(0,2*pi,100)'; 
circsx = p.*cos(t) ; 
circsy = p.*sin(t) ; 

if (fig1)
	plot(circsx,circsy);
end

%x = linspace(-1,1,100)';
%noise=normrnd(0,sigma,100,1);
%y = -x*cot(theta) + csc(theta) * [ p + noise ];
%y = -x/tan(theta) + 1/sin(theta) * [ p + noise ];

%if (fig1)
figure(1)
	plot(x,y, 'g.');
%end

% the slope and the intercept themselves do not come from a distribution
% they don't need priors, we need a prior for the noise only

% likelihood

c_i = 1:100;
c_j = 1:100;
p_r = (c_i-75)/500+p;
theta_r = (c_j-25)/500+theta;
n=length(y);	
for i = c_i;
for j = c_j;
	err_i = y * sin(theta_r(j)) + x * cos(theta_r(j)) - p_r(i);
%	ll(i,j)=(sigma^-n) * exp(-1/(2*sigma) * err_i'*err_i);
	lll(i,j)=-n*log(sigma) + (-(2*sigma)^(-2) * err_i'*err_i);
	ll(i,j)=exp(lll(i,j));
end
end

% normalize, or else I've to use gnuplot and can't rotate
mm=max(max(ll));
ll = ll/mm;

% somehow the maximum of the likelihood function does not respond exactly to (theta, p)
% it is consequently a bit lower, at around 4.9, it doesn't matter what theta is
% remarkably with p=-5 it is at around -5.1 so also shifted with 0.1 down
% this is even the case with 0, so probably I've made a mistake somewhere
figure(2)
surf(p_r,theta_r,ll(c_i,c_j))

