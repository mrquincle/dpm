stat=zeros(3,4);

plot_results='line';
plot_results='segment';

switch(plot_results)
case 'line'
    path='aug-30/';
    plot_title = 'Clustering performance of line detection algorithm';
    ofile = 'barplot-lines.png';
case 'segment'
    path='may-04/';
    plot_title = 'Clustering performance of segment detection algorithm';
    ofile = 'barplot-segments.png';
end

for i=1:3
	switch(i)
	case 1
		suc=importdata(strcat(path,'ri.avg.txt'));
	case 2
		suc=importdata(strcat(path,'ar.avg.txt'));
	case 4
		suc=importdata(strcat(path,'mi.avg.txt'));
	case 3
		suc=importdata(strcat(path,'hi.avg.txt'));    
	end

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

fg=figure(1);

errorbar(stat(:,1), stat(:,2), stat(:,3), stat(:,4), '#~.r')
ylim([0, 1]);
%set(gca,'XTickLabel',{'','Rand Index', 'Adjusted Rand', 'Mirkin', 'Hubert'})
set(gca,'XTickLabel',{'','Rand Index', '', 'Adjusted Rand', '', 'Hubert'})
title(plot_title);
%legend('Average performance');

W = 4; H = 3;
set(fg,'PaperUnits','inches');
set(fg,'PaperOrientation','portrait');
set(fg,'PaperSize',[H,W]);
set(fg,'PaperPosition',[0,0,W,H]);

FN = findall(fg,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(fg,'-property','FontSize');
set(FS,'FontSize',12);

saveas(1, ofile);


%print -dpng barplot.png
