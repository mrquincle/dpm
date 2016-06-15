mu=0;
gamma=[0.25 0.5 1 2];
t=-pi:0.01:pi;

for i=1:4
    S.gamma = gamma(i);
    S.mu = mu;
    y=cauchycdf(S,t);
    plot(t,y);
    axis([-pi pi 0 1]);
    hold on;
end

