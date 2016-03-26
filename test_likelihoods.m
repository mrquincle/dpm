
addpath('../inference/likelihood/normal');

method='DPM_Seg';

N=3;

mu0=[1 2];
mu1=[1 2.1];
mu2=[2 3];
mu=[mu0' mu1' mu2'];

% create x between a0 and b0
a0=5;
b0=10;
x=rand(1,1)*a0+b0-a0;
x=repmat(x, 1, N);
X=[ones(1,N);x];

mut=repmat(mu0, N, 1)';

y=sum(X.*mu,1);
data=[X; y];

b=[10 100 1000];

size(data)

clear R
N=1;
for i = 1:N
	R(i).mu = mu(:,1);
	R(i).Sigma = 1;
	R(i).a=0;
	R(i).b=b(i);
end

n=likelihoods(method, data, R)
