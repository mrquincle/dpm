

% define two Gaussian pdf's


% generate random variables from the two Gaussian pdf's

rv1=
rv2=

%%% this is not the difficult part, it is if rv1 and rv2 ends up mixed up in one vector
%%% in that case we don't know which random variable belongs to which component
%%% if we actually calculate the product, we don't know the form of the prior


% set hyper parameters for NIW product



% update NIW product with observations



% print new NIW parameters and compare with definition at start
