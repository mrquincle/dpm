clear all;
addpath("../utils");

D=2;

N = 3;

for i = 1:N
	S(i).M = rand(1,D);
	S(i).nu = rand(D,1);
	S(i).Sigma = rand(D,D);
end

printf "M is ordered as a vector\n"
M1 = S(1).M
tic
SM = extract(S, 'M');
toc
size(SM)

printf "nu is ordered as a column, extract goes much faster\n"
nu1 = S(1).nu
tic
Snu = extract(S, 'nu', 1);
toc

size(Snu)

printf "Sigma (DxD) can be extracted similarly"
%Sigma1 = S(1).Sigma
%Sigma2 = S(2).Sigma
%SSigma = extract(S, 'Sigma', 1);

%SSigma
%SM
%Snu
