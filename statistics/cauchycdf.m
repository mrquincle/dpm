% Just from Wolfram Alpha
% y=1/pi*atan(coth(gamma/2)*tan((t-mu)/2)) + 1/2;
function y = cauchycdf(S,t)
    y=1/pi*atan(coth(S.gamma/2)*tan((t-S.mu)/2)) + 1/2;
end

