% Generate samples from the Von Mises distribution.
%
% -- Function: x = vonmisesrnd(S)
%     Return a single sample from a Von Mises distribution
%     described by its sufficient statistics S. The sufficient statistics
%     consist out of:
%     * S.mu: location vector
%     * S.K >= 0: concentration parameter
% -- Function: x = vonmisesrnd(S, N)
%     Return N samples from the Von Mises distribution described by
%     its sufficient statistics S.
%
% For more information, see:
%  * https://en.wikipedia.org/wiki/Von_Mises_distribution
%
% The implementation for the generation of these random variates:
%  * Sample \sigma^2 from an inverse gamma distribution with parameters a and b
%  * Sample x from a normal distribution with mean \mu and variance \sigma^2/\lambda

function x = vonmisesrnd(S, N)
    x = vmrnd(S.MU, S.K, N);
%	if ~exist('N','var')
%		N=1;
%	end

    % for K approaching zero, sample uniform from [0, 2pi]
 %   if (S.K < 1e-6) 
  %      x = 2*pi*rand(N,1);
   %     return;
    %end
end
