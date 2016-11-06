clear all
close all

addpath('..');

display_in_2d = false;
display_cumulative = false;

a = 0; 
b = 24; 
 
mu1 = 6; 
kappa1 = 6; 
mu2 = 21; 
kappa2 = 8; 

weights0 = 0;
weights1 = 0.6;
weights2 = 1 - weights1 - weights0;

sum_weights = sum([weights0, weights1, weights2])

% difference between subsequent samples is 0.1
D=0.1;

C = [a:D:b];

if (display_in_2d)
	[C1 C2] = meshgrid(C);
	X = [C1(:) C2(:)];
else
	X = C;
end

N = size(X, 1);

Y = schedulepdf(X, weights0, weights1, weights2, a, b, mu1, kappa1, mu2, kappa2);

if (display_in_2d)
	X1 = X(:,1);
	X2 = X(:,2);

	yy = reshape(Y, size([C1]));

	% integrate over dimension 1
	y1 = cumsum(yy, 1) * D;
	y2 = cumsum(y1, 2) * D;
else
	yy = Y;
	y1 = cumsum(yy, 2) * D;
end

j=1;
f(j) = figure;

if (display_cumulative)
	if (display_in_2d) 
		mesh(C1, C2, y2);
	else
		area(X, y1);
		axis([0 24]);
	end
else
	if (display_in_2d) 
		mesh(C1, C2, yy);
	else
		%plot(X, Y);
		area(X, Y);
		axis([0 24]);
	end
end

% integrate this 

if (display_in_2d)
	S=sum(sum(yy)) * D * D;
else
	S=sum(yy) * D * pi / 12;
end

file=['image' num2str(j)];
saveas(f(j), file, 'png')

S

% S should be almost equal to 1
assert(S-1 < 0.01);

