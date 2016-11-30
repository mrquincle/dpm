% -- Function File: S = wishrnd (N, LAMBDA)
%     Draw S from a Wishart distribution with parameters NU and LAMBDA.
%     
%     The parameter N is the degrees of freedom (N > p - 1) 
%     The parameter LAMBDA > 0 is a scale matrix of dimension p x p
%
function W = wishrnd(n,lambda)

    [p p2]=size(lambda);

    if p~=p2
        warning('The matrix lambda to generate a Wishart distribution should be square');
    end

    if n <= (p - 1)
        warning('The number of degrees should equal the matrix dimension minus one or larger');
    end

    x = chol(lambda)' * randn(p,n);
    W=x*x';

end
