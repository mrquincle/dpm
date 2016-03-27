
% Update sufficient statistics of Normal Wishart distribution by removing
% data point z

mu0 = S.mu;
kappa0 = S.kappa;
nu0 = S.nu;
lambda0 = S.lambda;

kappa1=kappa0-1;
nu1=nu0-1;
mu1=kappa0/(kappa0-1)*mu0-1/(kappa0-1)*z;
lambda1=lambda0-(kappa0-1)/kappa0*(z-mu1)*(z-mu1)';
% if isnan(lambda1)
%     keyboard
% end
S.kappa = kappa1;
S.mu = mu1;
S.nu = nu1;
S.lambda = lambda1;
