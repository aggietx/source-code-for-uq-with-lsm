function [rmstime,rmsall,ccall ]=compute_velerrorGBsw(u_hat,gamma_mean,rk,kk,L1,Dim_Grid,st)
T=size(u_hat,2);
allux=zeros(Dim_Grid^2,T);alluy=zeros(Dim_Grid^2,T);
alluestx=zeros(Dim_Grid^2,T);alluesty=zeros(Dim_Grid^2,T);
rmstime=zeros(T,1);

for time=1:T
    [ux,uy]=computeGBvel(rk,kk,L1,Dim_Grid,u_hat(:,time));
    [uxest,uyest]=computeGBvel(rk,kk,L1,Dim_Grid,gamma_mean(:,time));
    allux(:,time)=ux(:);alluy(:,time)=uy(:);
alluestx(:,time)=uxest(:);alluesty(:,time)=uyest(:);
[rmsveltime]=computermsvel(ux,uy,uxest,uyest,0);
rmstime(time)=rmsveltime;
end


rmsall=zeros(Dim_Grid^2,1);
ccall=zeros(Dim_Grid^2,1);
for i=1:Dim_Grid^2
 [rmsx,ccx]=  rmscc(alluestx(i,st:end),allux(i,st:end),0);
  [rmsy,ccy]=   rmscc(alluesty(i,st:end),alluy(i,st:end),0);
  rmsall(i)=rmsx/2+rmsy/2;
  ccall(i)=ccx/2+ccy/2;
end
% % fprintf('Kmax is  %d\n',Kmax);

fprintf('mean (over time) vel  rms  %2.2f  \n',mean(rmstime(st:end)));
fprintf('mean (over point) vel normalized rms, cc are  %2.2f  %2.2f\n',mean(rmsall(1:end)),mean(ccall(1:end)));
fprintf('grid size is  %d\n',Dim_Grid);

