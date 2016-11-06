% Plots the probability density function of the Pareto distribution

addpath('../inference')
addpath('../inference/prior/pareto')

alpha=5;
N=1000;
lambda_r=0;
lambda_l=-4;

% this is not what we want, it would lead to a different shape for the
% pdf on the right side versus the left side
%p1=rnd_pareto(lambda_r, alpha, 1, N);
%p2=rnd_pareto(lambda_l, alpha, 1, N);
%p=[p1 p2];

[p1 p2]=rnd_pair_two_sided_pareto(lambda_l, lambda_r, alpha, 0, N);

p=[p1; p2];

fg=figure(1);

hist(p,1000);

i_dont_remember=false;

if (i_dont_remember)
    x=[0:0.001:5];

    k=3;
    b=1;

    % p(x) = k*b^k / x^(k+1) for x < b

    y=(x > b) * (k*b^k) ./ x.^(k+1);

    plot(x,y,'b-');

    hold on

    k=2;
    b=2;
    y=(x > b) * (k*b^k) ./ x.^(k+1);

    plot(x,y,'g-');
end

W = 7; H = 3;
set(fg,'PaperUnits','inches')
set(fg,'PaperOrientation','portrait');
set(fg,'PaperSize',[H,W])
set(fg,'PaperPosition',[0,0,W,H])

%ha = axes('Position', [0 0 1 1], 'Xlim', [0 1], 'Ylim', [0 1], 'Box', 'off', 'Visible', 'off', 'Units', 'normalized', 'clipping', 'off');
%text(0.5, 0.98,'Update of Pareto distribution after observations','HorizontalAlignment','center','VerticalAlignment', 'top')

FN = findall(fg,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(fg,'-property','FontSize');
set(FS,'FontSize',12);

saveas(1, "pareto_rnd.png");

