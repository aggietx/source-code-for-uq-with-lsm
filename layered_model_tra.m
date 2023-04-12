%%%% layered model
%%%% generate trajectory
clear;close all;rng(10);
L=24;
p=.5*2;rk=0;if p==1;testcase=1;else testcase=2;end
T=4000;dt=1/1000;N=T/dt;L1=2*pi;fprintf('L is %d \n',L);
x=zeros(L,N);sig_ex=.25;
y=x;gap=1;dt0=1*dt;N0=T/dt0;
x(:,1)=L1*rand(L,1);
y(:,1)=L1*rand(L,1);
Kmax=6;%theta=pi/3;
fprintf('Kmax is %2.2f\n',Kmax);
fprintf('p is %2.2f\n',p);
fprintf('sig ex is %2.2f\n',sig_ex);
kk0=(1:Kmax)';lx=1;
dk=0.0125*ones(Kmax,1);du=0.0125;H1=1;H2=1/2;
beta=2;
sigmak=zeros(Kmax,1);hk=zeros(Kmax,1);
sigmau=1/20/sqrt(2);sigmak(1:2)=1/20/sqrt(2);
hk(1)=H1/2-H1/2*1i;hk(2)=H2/2-H2/2*1i;

for ik=3:Kmax
    sigmak(ik)=1/20/sqrt(2)/(ik)^p;
    theta=rand(1)*2*pi;
    hk(ik)=sin(theta)/4/(ik)^p-cos(theta)/4/(ik)^p*1i;
    
end
kk1dfull=zeros(2,2*Kmax);
for i=1:Kmax
   kk1dfull(1,2*i-1)=i;
   kk1dfull(1,2*i)=-i;
end
kk1dfullf=[[0;0],kk1dfull];
repk1d=repmat(kk1dfull(1,:),L,1)*1i;
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
x_loc = [x(:,i-1),y(:,i-1)];     
if rk~=1
% %     psik(:,i)=psik(:,i-1)-dk.*psik(:,i-1)*dt+...
% %         1i*kk0*lx.*(beta./(kk0.^2)/lx^2-u(:,i-1)).*psik(:,i-1)*dt+...
% %         1i*kk0*lx./(kk0.^2)/lx^2.*hk.*u(:,i-1)*dt+sigmak.*randn(Kmax,1)*sqrt(dt);
    psik(:,i)=psik(:,i-1)-dk.*psik(:,i-1)*dt+...
        1i*kk0*lx.*(beta./(kk0.^2)/lx^2-u(:,i-1)).*psik(:,i-1)*dt+...
        1i*kk0*lx./(kk0.^2)/lx^2.*hk*u(:,i-1)*dt+1/sqrt(2)*sigmak.*(randn(Kmax,1)+1i*randn(Kmax,1))*sqrt(dt);    
    psistar=real(psik(:,i-1))-1i*imag(psik(:,i-1));
    u(:,i)=u(:,i-1)-du*u(:,i-1)*dt-1i*lx*sum(kk0.*hk.*psistar)*dt+...
        -1i*lx*sum(-kk0.*hkstar.*psik(:,i-1))*dt+sigmau*randn*sqrt(dt);
else
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
   
end
%% compute vel
a=psik(:,i-1);b=conj(a);%real(a)-1i*imag(a);
a1=1i.*a.*kk0;b1=-1i*b.*kk0;
uy1d = exp(1i * x(:,i-1) * kk0')*a1+exp(-1i * x(:,i-1) * kk0')*b1;
    if max(abs(imag(uy1d)))>10^(-6) 
      disp('complex velocity'  )
      max(abs(imag(u)))
    end
    uy1d=real(uy1d);
ux=u(i-1);
uy=uy1d+u(i-1);

% Gy1d = [zeros(L,1),exp(1i * x_loc *  kk1dfull).*repk1d ];
% Gu=[ones(L,1),zeros(L,2*Kmax)];
% G1=Gu;G2=(Gu+Gy1d);   
% fullunknow=zeros(1+2*Kmax,1);
% fullunknow(1)=u(i-1);fullunknow(2:2:end)=psik(:,i-1);
% fullunknow(3:2:end)=conj(psik(:,i-1));
% velx=G1*fullunknow;vely=G2*fullunknow;
% 
% if norm(ux-velx)>10^(-5) || norm(uy-vely)>10^(-5)
%     disp('error')
% end
%% update position
    x(:,i) = x(:,i-1) + ux * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in x
    y(:,i) = y(:,i-1) + uy * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in y
    x(:,i) = mod(x(:,i),L1);
    y(:,i) = mod(y(:,i),L1);    
    if max(abs(real(psik(:,i))))>1000
        disp('error,blow up');
        i
        return
    end
end
uhate=zeros(1+2*Kmax,T);
uhate(1,:)=u(:,1/dt:1/dt:end);
uhate(2:2:end,:)=psik(:,1/dt:1/dt:end);
uhate(3:2:end,:)=conj(psik(:,1/dt:1/dt:end));
xe=x(:,1/dt:1/dt:end);
ye=y(:,1/dt:1/dt:end);
xye=[xe;ye];
%  estuhat=zeros(2*Kmax+1,400);
% estuhat(1,:)=u(:,1/dt:1/dt:end);
%  estuhat(2,:)=psik(1,1/dt:1/dt:end); 
% estuhat(3,:)=conj(psik(1,1/dt:1/dt:end));
% estuhat(4,:)=psik(2,1/dt:1/dt:end);      
% estuhat(5,:)=conj(psik(2,1/dt:1/dt:end));
% K_max=10;uexact=u(:,1/dt*gap:1/dt*gap:end);
% psikexact=psik(:,1/dt*gap:1/dt*gap:end);
% compute_errorlayer1d
%  errest=[mean(rmstime(10:end)),mean(rmsall),mean(ccall)];
%  clear x y u psik G1 G1p G2 G2p Gu Gy1d

%  estuhat=zeros(2*Kmax+1,400);
% estuhat(1,:)=u(:,1/dt:1/dt:end);
% K_max=6;
% for i=1:K_max
%     estuhat(2*i,:)=psik(i,1/dt:1/dt:end); 
%      estuhat(2*i+1,:)=conj(psik(i,1/dt:1/dt:end)); 
% end
% K_max=10;uexact=u(:,1/dt*gap:1/dt*gap:end);
% psikexact=psik(:,1/dt*gap:1/dt*gap:end);
% compute_errorlayer1d
%   errest=[mean(rmstime(10:end)),mean(rmsall),mean(ccall)];
return
nx=30;ny=nx;Dim_Grid = nx;L1=2*pi;
x1d=linspace(0,L1,nx);x1d=x1d';
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
% return
close all;
for i=1/dt:5/dt:N; 
a=psik(:,i);b=real(a)-1i*imag(a);
a1=1i.*a.*kk0;b1=-1i*b.*kk0;
uy1d = exp(1i * x1d * kk0')*a1+exp(-1i * x1d * kk0')*b1;
uy=repmat(uy1d',ny,1)+u(i);
ux=ones(ny,nx)*u(i);
% 
    quiver(xx, yy, ux, uy, 'linewidth',1)
    hold on
    scatter(x(:,i),y(:,i),'red','linewidth',2);
    xlim([0, L1 ])
    ylim([0, L1 ])  
    xlabel('x','FontSize',18);  ylabel('y','FontSize',18)
    setgca(18)
    box on    
    title(['t = ', num2str(i*dt)],'FontSize',24)
    pause();
close all
end

return
