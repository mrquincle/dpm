% Plots the probability density function of the Pareto distribution

clear
clf

addpath('../inference/likelihood');

fig=figure(1);

x=[-5:0.001:20];

% x_m
k=[2 6];

% x_n
l=[5 9];

% alpha
b=[3 3];

% Number of plots
M=2;

% Populate matrix (multiple y-arrays)
y = zeros(M,length(x));

set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(fig)))

hold all
for i = 1:M
	%y(i,:)=paretopdf(x, b(i), k(i));
	y(i,:)=biparetopdf(x, b(i), k(i), l(i));

	plot(x,y(i,:),'-');
end

%print -dpng conflicting_segments.png
