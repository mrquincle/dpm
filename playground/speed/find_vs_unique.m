
% When M >> N, then unique is faster
% If N >> M (our case), then find is faster
N = 10000;

M = 2000000;
m = zeros(1, M);
c = zeros(N, 1);

for i=1:N
	k = ceil(M*rand);
	c(i) = k;
	m(k) = m(k) + 1;
end

tic
inda = unique(c);
toc


tic
indb = find(m);
toc

