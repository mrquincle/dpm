%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We generate a line for which we use theta and p as parameters. Here theta denotes the rotation
% of the line and p the distance of the point on the line that is closest to the origin. This is
% called the normal form of a line.
%
% The points on the line are distributed ortogonal to the line itself, as assumed in total least 
% squares, Deming regression, or errors-in-variables methods. It is not just in the y-direction as 
% is assumed in ordinary linear regression.
%
% We make the line into a segment, by "traveling" over the line starting from the point that is
% closest to the origin. Let this point by LP, then we sample points over the line 
% (orientation theta + 90 degrees) between d1 and d2. For example with theta=0, p=2, d1=4, d2=5, 
% there will be a vertical line at x=2 with points between y=4 and y=5. To have the same line
% segment at x=-2, the parameters are theta=pi, p=2, d1=-4, d2=-5. The distance p is always 
% nonnegative.
%
% For all possible values of theta and p (with a certain resolution)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear figures and variables
clf
clear

% Parameters that can be adjusted

% angle w.r.t. origin
theta=rand*2*pi;
% distance with origin
p=unifrnd(-10,10);
% noise
sigma=0.01;
% endpoints of the segment
d1=5; % travel d1 over line from LP, this is endpoint 1
d2=8; % travel d2 over line from LP, this is endpoint 2
% number of points
N=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualization parameters

% do we draw figure 1?
fig1=false;
% do we draw figure 2?
fig2=true;

% resolution of x and y coordinate (more than 100 takes a long time)
resolution=[100 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print parameters
printf("theta: %d\n", theta);
printf("p: %d\n", p);
printf("sigma: %d\n", sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate points using the angularrnd function
S.sigma = sigma;
T.theta = theta;
T.p = p;
T.d1 = d1;
T.d2 = d2;
[x y] = angularrnd(S, T, N);

if (fig1)
    figure(1)
	plot(x,y, 'g.');
end

% the slope and the intercept themselves do not come from a distribution
% they don't need priors, we need a prior for the noise only

% likelihood

% For different transformation parameters calculate the likelihood. So, we sweep over Tl and adjust
% Tl.theta and Tl.p.


c_i = 1:resolution(1);
c_j = 1:resolution(2);
p_arr = (c_i-resolution(1)/2)/50+p;
theta_arr = (c_j-resolution(2)/2)/50+theta;
for i = c_i;
for j = c_j;
    Tl.p = p_arr(i);
    Tl.theta = theta_arr(j);
%	err_i = y * sin(theta_r(j)) + x * cos(theta_r(j)) - p_r(i);
%	ll(i,j)=(sigma^-n) * exp(-1/(2*sigma) * err_i'*err_i);
%	lll(i,j)=-n*log(sigma) + (-(2*sigma)^(-2) * err_i'*err_i);
    lll(i,j)=logangularpdf(S, Tl, x, y);
	ll(i,j)=exp(lll(i,j));
end
end

% normalize, or else I've to use gnuplot and can't rotate
mm=max(max(ll));
lnorm = ll/mm;

% somehow the maximum of the likelihood function does not respond exactly to (theta, p)
% it is consequently a bit lower, at around 4.9, it doesn't matter what theta is
% remarkably with p=-5 it is at around -5.1 so also shifted with 0.1 down
% this is even the case with 0, so probably I've made a mistake somewhere
if (fig2)
    figure(2)
    surf(theta_arr,p_arr,lnorm(c_i,c_j))
    view(2)
end

