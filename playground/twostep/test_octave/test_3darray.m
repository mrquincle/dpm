
addpath("../../../utils")

K = 3;

dataset = {};
for i = 1:K
	dataset(i).data = [i i; i + 0.2 i + 0.2; i + 0.4 i + 0.4];
end

% Input
dataset.data

% Output
extract(dataset, 'data')

