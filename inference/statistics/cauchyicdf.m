% y=1/pi*atan(coth(gamma/2)*tan((t-mu)/2)) + 1/2;
% then calculate inverse for t
% useful for inverse transform sampling
function t = cauchyicdf(S,y)
    t = 2*atan(tanh(S.gamma/2)*tan(pi*(y-1/2)))+S.mu;
end


