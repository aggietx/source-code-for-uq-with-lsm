%  ix=1;iy=1;
%  i1=getind(ix,iy,kk);i2=getind(-ix,-iy,kk);
%  s=[i1,i2];
%   i1=getind(ix,iy,kkold);i2=getind(-ix,-iy,kkold);
%  s1=[i1,i2];

 
step=2;
% data=gamma_mean_trace(:,step);
% data=utest(:,50);
data=u_hat(:,step);
% data=sigmafit*randn(size(u_hat,1),1);
for ix=-Kmax:Kmax;
for iy=0:Kmax
 i1=getind(ix,iy,kk);i2=getind(-ix,-iy,kk);
 s=[i1,i2];
 if abs(imag((sum(data(s)  ))))>10^(-7) %|| abs(real(data(s(1)))-real(data(s(2))))>10^(-7)
     disp('error, not symmetric')
     data(s)
%      Ffit(s)
%      sigmak2(s)
     pause()   
     
 end
%  Ffit(s)
end
end

