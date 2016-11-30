clear all;
close all;

data_dir = '.';
data_dir = 'generated';

data_includes_labels=true;

data_glob=[data_dir '/*.data.txt' ];
fileList = glob(data_glob);

file = fileList{1,1};

data=load(file);

print_struct_array_contents(true);

if (data_includes_labels)
	dataR = data.schedule_gen(:,1:end-1);
	dataL = data.schedule_gen(:,end);
else
	dataR = data.regulars;
	labels = [4 5 1 3 3 2 4 4 1 1 5 1 1 2 1 2 1 3 1 5 2 1 1 4 2 3 3 4 2 1 2 1 3 4 3 1 2 1 1 2];
	dataL = labels';
end
	
M = [dataR dataL];

[values, order] = sort(M(:,end));

R = M(order,:);

T = R(:,1:end-1);

K = size(T,1);

for k=1:K
	L = values(k);
	f(L) = figure(L);

	show_schedule(T(k,:),k);
	hold on;
end

J = unique(values);
for j = 1:length(J);
	figure(j);
	file=[data_dir num2str(j)];
	saveas(f(j), file, 'png');
end

close all;
