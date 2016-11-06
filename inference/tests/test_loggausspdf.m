
addpath('..');
addpath('../inference/likelihood/normal');

mu=[0 0];
sigma = [0.5 0; 0 2];

% difference between subsequent samples is 0.1
D=0.2;

C = [-16:D:16];

[C1 C2] = meshgrid(C);
X = [C1(:) C2(:)];

N = size(X, 1);

MU = repmat(mu, N, 1);
%SIGMA = repmat(sigma, N, 1);
SIGMA = sigma;

Y = exp(loggausspdf(X, MU, SIGMA));

X1 = X(:,1);
X2 = X(:,2);

yy = reshape(Y, size([C1]));

% integrate over dimension 1
y1 = cumsum(yy, 1) * D;
y2 = cumsum(y1, 2) * D;

mesh(C1, C2, y2)

% integrate this 

S=sum(sum(yy)) * D * D;
S

% S should be almost equal to 1
assert(S-1 < 0.0001);
