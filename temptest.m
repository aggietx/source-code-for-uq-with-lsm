N=100000;
Af=zeros(1000,N);
Ab=zeros(1000,N);
for i=N:-1:1
    tmp=randn(1000,1)+1i*randn(1000,1);
    Af(:,i)=tmp;
    Ab(:,N-i+1)=tmp;
end