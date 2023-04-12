function s=ifft2my(a)

[M,N]=size(a);



s=zeros(M,N);

%%%% M,N are odd
for m=0:M-1
    for n=0:N-1
        temp=0;
        for k=0:M-1
            for l=0:N-1
                temp=temp +a(k+1,l+1)*exp(-sqrt(-1)*2*pi*(m*k/M+n*l/N));
            end
        end
        s(m+1,n+1)=temp;
    end
end
s=s/M/N;