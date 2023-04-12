% fft([2 4 -1 6])
%%% test fft in matrix
F=fft(eye(2*K_max+1));
fft2mat = kron(F,F);

[xx,yy] = meshgrid(0:2*K_max, 0:2*K_max);
xy=[xx(:),yy(:)];
kkk=xy';
ft=exp(-1i * xy * kkk*2*pi/(2*K_max+1));%%%fft matrix

d1=randn(2*K_max+1,2*K_max+1);

fft2d1=fft2(d1);
ft*d1(:)-fft2d1(:);


invft=exp(1i * xy * kkk*2*pi/(2*K_max+1))/(2*K_max+1)/(2*K_max+1);%%%% inverse fft matrix
norm(invft-inv(fft2mat));


% % data=u_hat(:,100);data1=data(1:(2*K_max+1)^2);
% % invG=exp(1i * kk0' * kk0/(2*K_max+1)*2*pi)
% % G=exp(-1i * kk0' * kk0/(2*K_max+1)*2*pi)
% % norm(G*invG/(2*K_max+1)^2-eye((2*K_max+1)^2))

data=u_hat(:,1000);
xy=kk0'/(2*K_max+1)*2*pi;
invG=exp(1i * xy*kk0);%%% inverse shallow water matrix
G=exp(-1i  * xy*kk0)/(2*K_max+1)^2;%%%%% shallow water matrix from physical domain to F domain
norm(G*invG-eye((2*K_max+1)^2));

invG1=invG .* (ones(length(xy),1) * rk(1,1:dimuhat0)); 
invG2=invG .* (ones(length(xy),1) * rk(2,1:dimuhat0)); 
invG3=invG .*(ones(length(xy),1) * rketa(1:dimuhat0) );
    u = (invG1* data(1:dimuhat0) ); 
    v = (invG2* data(1:dimuhat0)); 
    eta=(invG3* data(1:dimuhat0));
    %%% transform to Fourier domain
Gu=G*u;
Gv=G*v;
data1=zeros(dimuhat0,1);
for id=2:dimuhat0
    ka=kk0(1,id);kb=kk0(2,id);
    if rk(1,id)==0
      data1(id)=Gv(id)/rk(2,id);  
    else
       data1(id)=Gu(id)/rk(1,id);   
    end
end
norm(data1-data(1:dimuhat0))