N=20; % test, generate assignments to 9 possible tables
c=floor(unifrnd(0,4,1,N))*2+3
y=(1:N)*10

[m,~]=hist(c,1:9)

ind=find(m)

find(c==ind(1))
