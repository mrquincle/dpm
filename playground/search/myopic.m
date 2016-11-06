% Quickly figure out how far Normal distributions have to be separated from each other to become unreachable due to floating point limitations.

format long e
x=1:50;
y=normpdf(x);

y

format
nonzero=length(y(y>0))

format long e
x=10.^(1:1:20);
y=lognpdf(x);

y

