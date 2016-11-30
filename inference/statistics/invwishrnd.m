% -- Function File: S = invwishrnd (NU, LAMBDA)
%     Draw S from an inverse Wishart distribution with parameters NU and LAMBDA.
%     
%     The parameter NU is the degrees of freedom (NU > p - 1) 
%     The parameter LAMBDA is a scale matrix of dimension p x p
%
%     Internally, generating an inverse Wishart distribution is done by
%     generating a Wishart distribution and inverting the resulting matrix.
%
function S = invwishrnd(nu,lambda)
	iS=wishrnd(nu,lambda^-1);
	S=inv(iS);
end
