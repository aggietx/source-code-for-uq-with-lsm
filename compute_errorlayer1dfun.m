function [meanrmstime,meanrms,meancc]=compute_errorlayer1dfun(estuhat,uhate,Kmax,K_max)
st=3;
T=size(estuhat,2);nx=30;ny=nx;Dim_Grid = nx;L1=2*pi;
kk1dfull=zeros(2,2*Kmax);
for i=1:Kmax
   kk1dfull(1,2*i-1)=i;
   kk1dfull(1,2*i)=-i;
end
% ix=min(K_max,k_0);iy=min(K_max,k_0);ind=getind(ix,iy,kk);%% sample plot
rmstime=zeros(T,1);
% redind=getreducedind(kk,K_max,Kuse);
redind=1:size(estuhat,1);
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
x_loc=[xx(:),yy(:)];
repk1d2=repmat(kk1dfull(1,:),Dim_Grid^2,1)*1i;
Gy1d = [zeros(Dim_Grid^2,1),exp(1i * x_loc *  kk1dfull).*repk1d2 ];
Gu=[ones(Dim_Grid^2,1),zeros(Dim_Grid^2,2*Kmax)];
G1=Gu;G2=(Gu+Gy1d);
[G1p,G2p]=get1dG(K_max,Dim_Grid,L1);
allux=zeros(Dim_Grid^2,T);alluy=zeros(Dim_Grid^2,T);
alluestx=zeros(Dim_Grid^2,T);alluesty=zeros(Dim_Grid^2,T);
for time=1:T;
    i=time;
% fullunknow=zeros(1+2*Kmax,1);
% fullunknow(1)=uexact(i);fullunknow(2:2:end)=psikexact(:,i);
% fullunknow(3:2:end)=conj(psikexact(:,i));
fullunknow=uhate(:,i);
ux=G1*fullunknow;allux(:,time)=ux;
uy=G2*fullunknow;alluy(:,time)=uy;
uxest=G1p*estuhat(:,time);alluestx(:,time)=uxest;
uyest=G2p*estuhat(:,time);alluesty(:,time)=uyest;

[rmsveltime]=computermsvel(ux,uy,uxest,uyest,0);
rmstime(time)=rmsveltime;
% ccall(time)=ccvel;
end
fprintf('Kmax is  %d\n',Kmax);
fprintf('mean (over time) vel  rms  %2.2f  \n',mean(rmstime(st:end)));
meanrmstime=mean(rmstime(st:end));
rmsall=zeros(Dim_Grid^2,1);
ccall=zeros(Dim_Grid^2,1);
for i=1:Dim_Grid^2
 [rmsx,ccx]=  rmscc(alluestx(i,st:end),allux(i,st:end),0);
  [rmsy,ccy]=   rmscc(alluesty(i,st:end),alluy(i,st:end),0);
  rmsall(i)=rmsx/2+rmsy/2;
  ccall(i)=ccx/2+ccy/2;
end
fprintf('mean (over point) vel normalized rms, cc are  %2.2f  %2.2f\n',mean(rmsall),mean(ccall));
meanrms=mean(rmsall);
meancc=mean(ccall);
fprintf('grid size is  %d  \n',Dim_Grid);
clear alluestx alluesty allux alluy
return
load ccall
figure();
subplot(2,1,1);
title('uexact','FontSize',18)
setgca(18);
 plot(uexact,'*-','linewidth',1.5)
subplot(2,1,2);
plot(ccall,'*-','linewidth',1.5)
setgca(18);
title('cc','FontSize',18)