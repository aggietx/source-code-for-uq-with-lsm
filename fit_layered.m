%%%% layered model, another case.
%%%% with direction (lx,ly)
clear;close all;
T=50000;dt=1/1000;N=T/dt;L1=2*pi;rk=1;
% maxlag0=500000;maxlag=maxlag0/100;funfit=0;ts=1;fprintf('Lag is %d %d \n',maxlag0,maxlag);
testcase=2;%%%1:non gauss 2:strong non gauss
Kmax=6;p=2;%theta=pi/3;
fprintf('p is %2.2f\n',p);
fprintf('test case is %d\n',testcase);

kk0=(1:Kmax)';lx=4/5;ly=sqrt(1-lx^2);fprintf('(lx, ly) are (%2.2f,%2.2f)\n',lx,ly);
dk=0.0125*ones(Kmax,1);du=0.0125;H1=1;H2=1/2;
beta=2;sigmak=zeros(Kmax,1);hk=zeros(Kmax,1);
sig_ex=.25;
gap=100;dt0=1*dt;N0=T/dt0;
fprintf('Kmax is %d\n',Kmax);
fprintf('sig ex is %2.2f\n',sig_ex);

if testcase==1;
    sigmau=1/2/sqrt(2);sigmak(1)=1/1/sqrt(2);sigmak(2)=1/2/sqrt(2);
else
    sigmau=1/1/sqrt(2);sigmak(1)=1/2/sqrt(2)/2;sigmak(2)=1/2/sqrt(2)/2;
end
hk(1)=H1/2-H1/2*1i;hk(2)=H2/2-H2/2*1i;
rng(10);
% a=real(psik);b=imag(psik);c=sum(a.^2+b.^2,2);c(1)/c(2)
for ik=3:Kmax
    if testcase==1
    sigmak(ik)=1/ik^p/sqrt(2);
    else
     sigmak(ik)=1/ik^p/sqrt(2)/2;   
    end
    theta=rand(1)*2*pi;
    hk(ik)=sin(theta)/4/(ik)^p-cos(theta)/4/(ik)^p*1i;
    
end
hkstar=real(hk)-1i*imag(hk);
psik=zeros(Kmax,N);u=zeros(1,N);
tempu=u(:,1);
temppsik=psik(:,1);

kk1dfull=zeros(2,2*Kmax);
for i=1:Kmax
   kk1dfull(1,2*i-1)=i*lx;
   kk1dfull(1,2*i)=-i*lx;
   kk1dfull(2,2*i-1)=i*ly;
   kk1dfull(2,2*i)=-i*ly;   
end

kk1dfullf=[[0;0],kk1dfull];
kx=kk1dfullf(1,:);ky=kk1dfullf(2,:);
dX = 1i*kx*2*pi/L1;
dY = 1i*ky*2*pi/L1;
dX=dX(:);dY=dY(:);

rng(10);
% sigmak=zeros(size(sigmak));sigmau=0;
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end
        
        %% step1
    k1a=-dk.*temppsik+1i*kk0*lx.*(beta./(kk0.^2)-tempu).*temppsik+...
        1i*kk0*lx./(kk0.^2).*hk.*tempu;
    k1b=-du*tempu-1i*lx*sum(kk0.*hk.* conj(temppsik))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik);
    %% step2
    temppsik2=temppsik+k1a*dt/2; tempu2=tempu+k1b*dt/2;
    k2a=-dk.*temppsik2+1i*kk0*lx.*(beta./(kk0.^2)-tempu2).*temppsik2+...
        1i*kk0*lx./(kk0.^2).*hk.*tempu2;
    k2b=-du*tempu2-1i*lx*sum(kk0.*hk.* conj(temppsik2))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik2);    
    %% step3
    temppsik3=temppsik2+k2a*dt/2; tempu3=tempu2+k2b*dt/2;
    k3a=-dk.*temppsik3+1i*kk0*lx.*(beta./(kk0.^2)-tempu3).*temppsik3+...
        1i*kk0*lx./(kk0.^2).*hk.*tempu3;
    k3b=-du*tempu3-1i*lx*sum(kk0.*hk.*conj(temppsik3))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik3); 
    %% step4
     temppsik4=temppsik3+k3a*dt; tempu4=tempu3+k3b*dt;
     k4a=-dk.*temppsik4+1i*kk0*lx.*(beta./(kk0.^2)-tempu4).*temppsik4+...
        1i*kk0*lx./(kk0.^2).*hk.*tempu4;
     k4b=-du*tempu4-1i*lx*sum(kk0.*hk.* conj(temppsik4))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik4); 
    ksa=1/6*k1a+1/3*k2a+1/3*k3a+1/6*k4a;
    ksb=1/6*k1b+1/3*k2b+1/3*k3b+1/6*k4b;
psik(:,i)=temppsik+ksa*dt+1/sqrt(2)*sigmak.*(randn(Kmax,1)+1i*randn(Kmax,1))*sqrt(dt);  
% % % psik(:,i)=temppsik+ksa*dt+(1+1i)/sqrt(2)*sigmak.*randn(Kmax,1)*sqrt(dt); 
% % psik(:,i)=temppsik+ksa*dt+sigmak.*(randn(Kmax,1))*sqrt(dt); 
u(:,i)=tempu+ksb*dt+sigmau*randn*sqrt(dt);
temppsik=psik(:,i);
tempu=u(:,i);
end
traj=[u;psik];

traj1=zeros(1+2*Kmax,N);traj1(1,:)=u;
traj1(2:2:end,:)=psik;traj1(3:2:end,:)=conj(psik);
kk1dfullfake=zeros(2,2*Kmax);
for i=1:Kmax
   kk1dfullfake(1,2*i-1)=i;
   kk1dfullfake(1,2*i)=-i; 
end
kk1dfullfake=[[0;0],kk1dfullfake];
