stat=zeros(4,4);

path='aug-30/';
%path='may-04/';

plot_title = 'Clustering performance of line detection algorithm'

for i=1:4
	switch(i)
	case 1
		suc=importdata(strcat(path,'ri.avg.txt'));
	case 2
		suc=importdata(strcat(path,'ar.avg.txt'));
	case 3
		suc=importdata(strcat(path,'mi.avg.txt'));
	case 4
		suc=importdata(strcat(path,'hi.avg.txt'));
	end

    suc

	stat(i,1)=i;
	avgsuc=mean(suc);
	stat(i,2)=avgsuc;
	ssuc=suc(suc<avgsuc);
	N=length(ssuc);
	serrsqsum=sum((ssuc-avgsuc).^2);
	serr=sqrt(serrsqsum/(N-1));
	stat(i,3)=serr;
	gsuc=suc(suc>avgsuc);
	N=length(gsuc);
	gerrsqsum=sum((gsuc-avgsuc).^2);
	gerr=sqrt(gerrsqsum/(N-1));
	stat(i,4)=gerr;
end
stat

errorbar(stat(:,1), stat(:,2), stat(:,3), stat(:,4), '#~.r')
ylim([0, 1]);
set(gca,'XTickLabel',{'','Rand Index', 'Adjusted Rand', 'Mirkin', 'Hubert'})
title(plot_title);
legend('Average performance');

print -dpng barplot.png
