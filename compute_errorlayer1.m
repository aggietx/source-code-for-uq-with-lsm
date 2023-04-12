%%% compute exact layer vel and estimated vel
%%% (lx,ly)
% get psikexact.mat  uexact.mat estuhat.mat 
%load psikexact; load uexact;load estuhat.mat 
% load data
st=10;
T=size(estuhat,2);nx=30;ny=nx;Dim_Grid = nx;L1=2*pi;
K_max=(sqrt(size(estuhat,1))-1)/2;Kuse=K_max;if K_max~=Kuse;reduced=1;else reduced=0;end;
kk0=(1:Kmax)';%dt=1/1000;%%%% for computing exact solution.
if Kuse>K_max;Kuse=K_max;reduced=0;end
dimuhat=(2*K_max+1)^2;
[ky,kx]=meshgrid([0:K_max,-K_max:-1],[0:K_max,-K_max:-1]);
kk=[kx(:),ky(:)]';
rk=zeros(size(kk));
for i=2:dimuhat
        ki=kk(:,i);
%         rk(:,i)=1i*[-ki(2);ki(1)]/norm(ki);
        rk(:,i)=1i*[-ki(2);ki(1)]/sqrt(ki(1)^2+ki(2)^2+1);
end
rk(:,1)=[1,1];
% ix=min(K_max,k_0);iy=min(K_max,k_0);ind=getind(ix,iy,kk);%% sample plot

% redind=getreducedind(kk,K_max,Kuse);
redind=1:dimuhat;
allux=zeros(Dim_Grid^2,T);alluy=zeros(Dim_Grid^2,T);
alluestx=zeros(Dim_Grid^2,T);alluesty=zeros(Dim_Grid^2,T);
rmstime=zeros(T,1);
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
    Ga = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(1,:))); % Fourier bases for u
    Gb = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(2,:))); % Fourier bases for v

for time=1:T;
    i=time;
x1d=linspace(0,L1,nx);x1d=x1d';
    uxest = (Ga* estuhat(:,time)); 
    uyest = (Gb* estuhat(:,time)); 
    if max(abs(imag(uxest)))>10^(-6) || max(abs(imag(uyest)))>10^(-6)
      disp('complex velocity, record exact'  )
      max(abs(imag(uxest)))
    end
    
    uxest=real(uxest);uyest=real(uyest);
 uxest=reshape(uxest, Dim_Grid,Dim_Grid);
 uyest=reshape(uyest, Dim_Grid,Dim_Grid);

temppsi=zeros(1+2*Kmax,1);
% temppsi(1,:)=uexact(i);
temppsi(2:2:end,:)=psikexact(:,i);
temppsi(3:2:end,:)=conj(psikexact(:,i));

ux=G1fixed*temppsi+uexact(i);
uy=G2fixed*temppsi+uexact(i);
if max(abs(imag(ux(:))))>10^(-6) || max(abs(imag(uy(:))))>10^(-6)
    disp('error, complex velicity')
end
ux=real(ux);uy=real(uy);
ux=reshape(ux, Dim_Grid,Dim_Grid);
uy=reshape(uy, Dim_Grid,Dim_Grid);

allux(:,time)=ux(:);alluy(:,time)=uy(:);
alluestx(:,time)=uxest(:);alluesty(:,time)=uyest(:);
% % [rmsvel,ccvel]=rmscc(sqrt(uxest.^2+uyest.^2),sqrt(ux.^2+uy.^2),0);
% % rmsall(time)=rmsvel;
[rmsveltime]=computermsvel(ux,uy,uxest,uyest,0);
rmstime(time)=rmsveltime;
end
fprintf('Kmax is  %d\n',K_max);

fprintf('mean (over time) vel  rms  %2.2f  \n',mean(rmstime(st:end)));

rmsall=zeros(Dim_Grid^2,1);
ccall=zeros(Dim_Grid^2,1);
for i=1:Dim_Grid^2
 [rmsx,ccx]=  rmscc(alluestx(i,st:end),allux(i,st:end),0);
  [rmsy,ccy]=   rmscc(alluesty(i,st:end),alluy(i,st:end),0);
  rmsall(i)=rmsx/2+rmsy/2;
  ccall(i)=ccx/2+ccy/2;
end
fprintf('mean (over point) vel normalized rms, cc are  %2.2f  %2.2f\n',mean(rmsall(1:end)),mean(ccall(1:end)));
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