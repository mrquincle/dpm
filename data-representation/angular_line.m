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

addpath("../inference/utils")

% Parameters that can be adjusted

% angle w.r.t. origin
theta=rand*2*pi;
% distance with origin
p=unifrnd(-10,10);

% noise
sigma=0.1;
sigma=0.001;
% endpoints of the segment
d1=5; % travel d1 over line from LP, this is endpoint 1
d2=8; % travel d2 over line from LP, this is endpoint 2
% number of points
N=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualization parameters

% do we draw figure 1?
fig1=true;
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
[x y, E0, E1] = angularrnd(S, T, N);

if (fig1)
    figure(1);
    plot(x,y, 'go');
    hold on;
    E=[E0; E1];
    plot(E(:,1),E(:,2),'r-');
    axis("square");
    axis("equal");
    hold off;
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

% This is not incorrect. Check two points. If we draw a line through it with different origin-distances and angles
% we will find a nice periodic pattern. This is because even the origin-distance increases, the angle increases as 
% well, compensating for the overall line fit. 
%
% If we want to have a single peak we can bring down the noise, with sigma=0.001 you will have a single peak.
% You can also increase the number of points, although check if values are not clipping to infinity.
%
if (fig2)
    figure(2)
    surf(theta_arr,p_arr,lnorm(c_i,c_j))
    view(2)
    save "../output/angular_line.txt"
    C = meshtable(theta_arr, p_arr, lnorm);
    save "../output/angular_line.lnorm.txt" C
    D = [theta_arr' lnorm];
    D = [0 p_arr; D];
%    D = [[1:length(theta_arr)]' lnorm];
%    D = [0 1:length(p_arr); D];
    csvwrite("../output/angular_line.plotly.csv", D)
%    csvwrite("../output/angular_line.plotly.xaxis.csv", theta_arr)
%    csvwrite("../output/angular_line.plotly.yaxis.csv", p_arr)
end

% get two largest values
[val,ind] = maximum(ll, 2);

[y,x] = inds2sub(size(ll),ind);

theta0 = theta_arr(x(1));
p0 = p_arr(y(1));

printf("Maximum theta found at: %d\n", theta0);
printf("Maximum p found at: %d\n", p0);

printf("Difference in theta found at: %d\n", theta0 - theta);
printf("Difference in p found at: %d\n", p0 - p);

theta1 = theta_arr(x(2));
p1 = p_arr(y(2));

printf("Second best theta found at: %d\n", theta1);
printf("Second best p found at: %d\n", p1);

printf("Difference in theta found at: %d\n", theta1 - theta);
printf("Difference in p found at: %d\n", p1 - p);

T
