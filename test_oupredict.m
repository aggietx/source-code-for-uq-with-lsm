%%%% conditioner gaussian test with fitted data
clear;
% load ../../../Downloads/allpsi1;
% close all
rng(10);
% load x1;load y1;
% load allpsi_k6_1year;allpsi1=allpsi_k6_1year;
% load x1year;load y1year;x1=x1year;y1=y1year;
% load dk;load omegak;load mu ;load sigmak2;

load allpsi1year;allpsi1=allpsi1year;%%% k=10;
load traj10;L=size(traj10,1)/4;x1=traj10(0*L+1:1*L,:);y1=traj10(1*L+1:2*L,:);
load dk_k10.mat ;load omegak_k10.mat;load mu_k10.mat ;load sigmak2_k10.mat
fprintf('L is %d\n',L);
sig_ex=0.02;Kmax=10;refine=1;50;N=1.5e4*refine;dt =1/refine* 2.5*1E-4*93*24*3600;%% seconds
L1=2*pi*400000;Nf=256;T0 = 93*24*3600;Len0 = 400000; % dimensionalized length scale
K_max=Kmax;fprintf('total time step is %d\n',N);fprintf('dt is %2.2f\n',dt);
dimuhat=(2*K_max+1)^2;dimuhatf=Nf^2;
[ky,kx]=meshgrid([0:K_max,-K_max:-1],[0:K_max,-K_max:-1]);
kk=[kx(:),ky(:)]';
kf = [0:Nf/2 -Nf/2+1:-1]';  % wavenumbers
[kyf,kxf]=meshgrid(kf,kf);
kkf=[kxf(:),kyf(:)]';

idx=zeros(2*Kmax+1,1);ind=1;
for i=-Kmax:Kmax
    idx(ind)=find(kf==i);ind=ind+1;
end
idy=idx;
kredx=kf(idx);kredy=kf(idy);
[kyr,kxr]=meshgrid(kredy,kredx);
kkr=[kxr(:),kyr(:)]';
        dX = 1i*repmat(kf',[Nf 1 2])*2*pi/L1;
    dY = 1i*repmat(kf,[1 Nf 2])*2*pi/L1;
    dX=dX(:,:,1);dY=dY(:,:,1);
dX=transpose(dX(:));
dY=transpose(dY(:));
redind1=[];
for iy=[0:K_max,-K_max:-1]
    for ix=[0:K_max,-K_max:-1]
        redind1=[redind1;getind(ix,iy,kkf)];
    end
end
dXr=dX(redind1);dYr=dY(redind1);
kkfr=kkf(:,redind1);
x=zeros(L,N);
y=zeros(L,N);
 x(:,1)=x1(:,1);
 y(:,1)=y1(:,1);
velmax=zeros(N,1);
 u_hat=zeros(dimuhat,N);
 sigmak=form_sigmatrix1((2*Kmax+1)^2,sigmak2,kk,Kmax);


for iy=-K_max:K_max 
    for ix=-K_max:K_max

ind=getind(iy,ix,kk);

u_hat(ind,1)=allpsi1(iy+K_max+1,ix+K_max+1,1);
% temp=(allpsi2(getind(iy,ix,kkr),1));norm(temp-u_hat(ind,1))
    end
end
iy=3;ix=3;ind=getind(iy,ix,kk);
a=u_hat(ind,:);
b=allpsi1(iy+K_max+1,ix+K_max+1,:);
plot2line(real(b(:)),real(a(:)));
title([' mode ( ', num2str(iy),' , ', num2str(ix), ' )'],'FontSize',16)

disp('generate signal...');
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end    
         x_loc = [x(:,i-1),y(:,i-1)];
            u_hat(:,i)=u_hat(:,i-1)+(-dk+1i*omegak).*u_hat(:,i-1)*dt+...
       mu*dt+sqrt(dt)*sigmak*(randn(dimuhat,1));
    G1 = -(exp(1i * x_loc * kkfr*2*pi/L1))/(Nf^2).* (ones(L,1) * dYr)*Len0/T0;
    G2 = (exp(1i * x_loc * kkfr*2*pi/L1))/(Nf^2).* (ones(L,1) * dXr)*Len0/T0;
   u=G1*u_hat(:,i);v=G2*u_hat(:,i);
if max(abs(real(u)))>10
    disp('blow up')
     return
end
   if max(abs(imag(u)))>10^(-5) || max(abs(imag(v)))>10^(-5) 
       disp('complex velocity');
       return
   else
       u=real(u);v=real(v);
   end
   velmax(i)=max(sqrt(u.^2+v.^2));
       x(:,i)=x(:,i-1)+u*dt+sqrt(dt)*sig_ex*randn(L,1);
    y(:,i)=y(:,i-1)+v*dt+sqrt(dt)*sig_ex*randn(L,1);

   
    x(:,i)=mod(x(:,i),L1);
    y(:,i)=mod(y(:,i),L1);

end

 return
% i=1;gap=10;
% scatter(x1(i,1:gap:end),y1(i,1:gap:end));hold on
% scatter(x(i,1:gap:end),y(i,1:gap:end));xlim([0,L1]);ylim([0,L1]);
outindx=[];
for i=1:L
    i
for j=1:size(x,2)-1
ind=0;
    if abs(x(i,j+1)-x(i,j))>L1/2 || abs(y(i,j+1)-y(i,j))>L1/2 || abs(x1(i,j+1)-x1(i,j))>L1/2 || abs(y1(i,j+1)-y1(i,j))>L1/2
    ind=ind+1;   
    end
 if ind>0
     break
 end
end   
if ind>0
 outindx=[outindx;i];
end
end
inindex=setdiff(1:L,unique(outindx));

figure()
gap=10;
for i=inindex(1:30);
plot(x1(i,1:gap:end)/1000,y1(i,1:gap:end)/1000,'k','linewidth',2);hold on
plot(x1(i,1)/1000,y1(i,1)/1000,'ob','linewidth',3);hold on
plot(x(i,1:gap:end)/1000,y(i,1:gap:end)/1000,'r','linewidth',2);xlim([0,L1]/1000);ylim([0,L1]/1000);
end
setgca(18)
return
% load traj;%%% k=6;
load traj10;traj=traj10;%% k=10;
gap=10;L=size(traj,1)/4;
x=traj(1:L,:);y=traj(1+L:2*L,:);%%% full
x1=traj(1+2*L:3*L,:);y1=traj(1+3*L:4*L,:);%%%% reduced k

outindx=[];
for i=1:L
for j=1:size(x,2)-1

    if abs(x(i,j+1)-x(i,j))>L1/2 || abs(y(i,j+1)-y(i,j))>L1/2 || abs(x1(i,j+1)-x1(i,j))>L1/2 || abs(y1(i,j+1)-y1(i,j))>L1/2
        outindx=[outindx;i];
    end
    
end   
end
inindex=setdiff(1:L,unique(outindx));

% close all;
figure()
for i=inindex(1:30);
plot(x1(i,1:gap:end)/1000,y1(i,1:gap:end)/1000,'r','linewidth',2);hold on
plot(x1(i,1)/1000,y1(i,1)/1000,'ob','linewidth',3);hold on
plot(x(i,1:gap:end)/1000,y(i,1:gap:end)/1000,'k','linewidth',2);xlim([0,L1/1000]);ylim([0,L1/1000]);hold on
end
setgca(18)
