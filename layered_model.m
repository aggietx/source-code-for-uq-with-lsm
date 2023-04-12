%%%% layered model
clear;close all;
T=20000;dt=1/1000;N=T/dt;
testcase=1;%%%1:non gauss 2:strong non gauss
Kmax=10;p=2;%theta=pi/3;
fprintf('p is %2.2f\n',p);
fprintf('test case is %d\n',testcase);
kk0=(1:Kmax)';lx=1;
dk=0.0125*ones(Kmax,1);du=0.0125;H1=1;H2=1/2;
beta=2;sigmak=zeros(Kmax,1);hk=zeros(Kmax,1);

if testcase==1;
    sigmau=1/2/sqrt(2);sigmak(1)=1/1/sqrt(2);sigmak(2)=1/2/sqrt(2);
else
    sigmau=1/1/sqrt(2);sigmak(1)=1/1/sqrt(2)/20;sigmak(2)=1/2/sqrt(2)/20;
end
hk(1)=H1/2-H1/2*1i;hk(2)=H2/2-H2/2*1i;
rng(10);
for ik=3:Kmax
    if testcase==1
    sigmak(ik)=1/ik/sqrt(2);
    else
     sigmak(ik)=1/ik/sqrt(2)/20;   
    end
    theta=rand(1)*2*pi;
    hk(ik)=sin(theta)/4/(ik)^p-cos(theta)/4/(ik)^p*1i;
    
end
hkstar=real(hk)-1i*imag(hk);
psik=zeros(Kmax,N);u=zeros(1,N);
tempu=u(:,1);
temppsik=psik(:,1);

rng(10);
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end
        %% step1
    k1a=-dk.*temppsik+1i*kk0*lx.*(beta./(kk0.^2)/lx^2-tempu).*temppsik+...
        1i*kk0*lx./(kk0.^2)/lx^2.*hk.*tempu;
    k1b=-du*tempu-1i*lx*sum(kk0.*hk.* conj(temppsik))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik);
    %% step2
    temppsik2=temppsik+k1a*dt/2; tempu2=tempu+k1b*dt/2;
    k2a=-dk.*temppsik2+1i*kk0*lx.*(beta./(kk0.^2)/lx^2-tempu2).*temppsik2+...
        1i*kk0*lx./(kk0.^2)/lx^2.*hk.*tempu2;
    k2b=-du*tempu2-1i*lx*sum(kk0.*hk.* conj(temppsik2))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik2);    
    %% step3
    temppsik3=temppsik2+k2a*dt/2; tempu3=tempu2+k2b*dt/2;
    k3a=-dk.*temppsik3+1i*kk0*lx.*(beta./(kk0.^2)/lx^2-tempu3).*temppsik3+...
        1i*kk0*lx./(kk0.^2)/lx^2.*hk.*tempu3;
    k3b=-du*tempu3-1i*lx*sum(kk0.*hk.*conj(temppsik3))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik3); 
    %% step4
     temppsik4=temppsik3+k3a*dt; tempu4=tempu3+k3b*dt;
     k4a=-dk.*temppsik4+1i*kk0*lx.*(beta./(kk0.^2)/lx^2-tempu4).*temppsik4+...
        1i*kk0*lx./(kk0.^2)/lx^2.*hk.*tempu4;
     k4b=-du*tempu4-1i*lx*sum(kk0.*hk.* conj(temppsik4))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik4); 
    ksa=1/6*k1a+1/3*k2a+1/3*k3a+1/6*k4a;
    ksb=1/6*k1b+1/3*k2b+1/3*k3b+1/6*k4b;
psik(:,i)=temppsik+ksa*dt+sigmak.*randn(Kmax,1)*sqrt(dt);
u(:,i)=tempu+ksb*dt+sigmau*randn*sqrt(dt);
temppsik=psik(:,i);
tempu=u(:,i);

%     psik(:,i)=psik(:,i-1)-dk.*psik(:,i-1)*dt+...
%         1i*kk0*lx.*(beta./(kk0.^2)/lx^2-u(:,i-1)).*psik(:,i-1)*dt+...
%         1i*kk0*lx./(kk0.^2)/lx^2.*hk.*u(:,i-1)*dt+sigmak.*randn(Kmax,1)*sqrt(dt);
%     psistar=conj(psik(:,i-1));
%     u(:,i)=u(:,i-1)-du*u(:,i-1)*dt-1i*lx*sum(kk0.*hk.*psistar)*dt+...
%         -1i*lx*sum(-kk0.*hkstar.*psik(:,i-1))*dt+sigmau*randn*sqrt(dt);
    
    if max(abs(real(psik(:,i))))>1000
        disp('error,blow up');
        i
        return
    end
end


return
acfpsi=zeros(Kmax+1,60/dt+1);
gap=20;psi=zeros(Kmax+1,1000/dt/gap);
for i=1:Kmax
    acfpsi(i,:) = autocorr(real(psik(i,:)),60/dt);
    temp=psik(i,1:1000/dt);
    psi(i,:)=temp(1:gap:end);
end
    temp=u(1:1000/dt);
    psi(end,:)=temp(1:gap:end);
acfpsi(end,:) = autocorr(u,60/dt);

nx=100;
streamfun=zeros(nx+1,1000);
x=0:2*pi/nx:2*pi;x=x';
for i=1:1000
a=psik(:,i/dt);b=real(a)-1i*imag(a);
streamfun(:,i) = exp(1i * x * kk0')*a+exp(-1i * x * kk0')*b;
end
streamfun75=streamfun;save streamfun75 streamfun75
 acfpsi1=acfpsi;save acfpsi1 acfpsi1;  
 psi1=psi;save psi1 psi1;
 
 
 data=u(1,:);
mu=mean(data);sig=var(data);
[~,r1]=ksdensity(data);
maxdata=max(r1);mindata=min(r1);
dx=(maxdata-mindata)/500;
xx=mindata:dx:maxdata;
[prob,r1]=ksdensity(data,xx);


datan=mu+sqrt(sig)*randn(2000000,1);
% % % [probn,rn]=ksdensity(datan);
maxdatan=max(datan);mindatan=min(datan);
dx=(maxdatan-mindatan)/500;
xxn=mindatan:dx:maxdatan;
  [probn,rn]=ksdensity(datan,xxn);
  acfv = autocorr(data,60/dt);
  
  data=[r1;prob;rn;probn];acfv=acfv;save data data; save acfv acfv

   figure();plot(data(1,:),data(2,:));hold on;plot(data(3,:),data(4,:),'--k')
 figure();plot(data(1,:),log(data(2,:)));hold on;plot(data(3,:),log(data(4,:)),'--k')
 figure();plot(acfv)
get data.mat acfv.mat



% % L1=2*pi;
% % Dim_Grid = 40;
% % [xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
% % xy=[xx(:),yy(:)];kk=[0;0];
% %     Ga = (exp(1i * xy * kk*2*pi/L1)); % Fourier bases for u
% %     Gb = (exp(1i * xy * kk*2*pi/L1)); % Fourier bases for v
% % u=Ga*0.0574;v=Gb*0.0574;
% % u=reshape(u,Dim_Grid,Dim_Grid);
% % v=reshape(v,Dim_Grid,Dim_Grid);
% % 
% % quiver(u,v)