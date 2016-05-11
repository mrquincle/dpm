stat=zeros(4,4);

for i=1:4
	switch(i)
	case 1
		suc=importdata('output/ri.avg.txt');
	case 2
		suc=importdata('output/ar.avg.txt');
	case 3
		suc=importdata('output/mi.avg.txt');
	case 4
		suc=importdata('output/hi.avg.txt');
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
stat

errorbar(stat(:,1), stat(:,2), stat(:,3), stat(:,4), '#~.r')
ylim([0, 1]);
set(gca,'XTickLabel',{'','Rand Index', 'Adjusted Rand', 'Mirkin', 'Hubert'})
title("Clustering performance of Algorithm 1 using several metrics")
legend("Average performance");
