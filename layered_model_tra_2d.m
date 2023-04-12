%%%% layered model
%%%% generate trajectory
clear;close all;rng(10);
L=24;
p=.5*2;rk=0;if p==1;testcase=1;else testcase=2;end
T=1000;dt=1/1000;N=T/dt;L1=2*pi;fprintf('L is %d \n',L);
x=zeros(L,N);sig_ex=.25;
y=x;gap=100;ts=1000;dt0=1*dt;N0=T/dt0;
x(:,1)=L1*rand(L,1);
y(:,1)=L1*rand(L,1);
Kmax=3;%theta=pi/3;
fprintf('Kmax is %2.2f\n',Kmax);
fprintf('p is %2.2f\n',p);
fprintf('sig ex is %2.2f\n',sig_ex);
dimuhat=(2*Kmax+1)^2;
[kx,ky]=meshgrid([0:Kmax,-Kmax:-1],[0:Kmax,-Kmax:-1]);
kx=kx(:);ky=ky(:);
kk=[kx(:),ky(:)]';
psik=zeros(dimuhat,N);
sigmak0=0.1*ones(dimuhat,1);sigmak0(1)=0;
sigmak=form_sigmatrix1(dimuhat,sigmak0,kk,Kmax);
sigmau=.1;
u=zeros(1,N);
beta=2;
dk=0.0125*ones(dimuhat,1);dk(1)=0;
du=0.0125;
ksquare=kx.^2+ky.^2;
% hk(1)=H1/2-H1/2*1i;hk(2)=H2/2-H2/2*1i;
% hk=zeros(Kmax,1);
% for ik=3:Kmax
%     sigmak(ik)=1/20/sqrt(2)/(ik)^p;
%     theta=rand(1)*2*pi;
%     hk(ik)=sin(theta)/4/(ik)^p-cos(theta)/4/(ik)^p*1i;
%     
% end
hk=ones(dimuhat,1);

k=[0:Kmax,-Kmax:-1]';
dX = 1i*kx*2*pi/L1;
dY = 1i*ky*2*pi/L1;
dX=dX(:);dY=dY(:);
rng(10);
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end
x_loc = [x(:,i-1),y(:,i-1)];     
if rk~=1
    psik(:,i)=psik(:,i-1)-dk.*psik(:,i-1)*dt+...
        1i*kx.*(beta./ksquare-u(:,i-1)).*psik(:,i-1)*dt+...
        1i*kx./ksquare.*hk*u(:,i-1)*dt+sigmak*randn(dimuhat,1)*sqrt(dt);    
    psistar=real(psik(:,i-1))-1i*imag(psik(:,i-1));
    u(:,i)=u(:,i-1)-du*u(:,i-1)*dt+1i*dt*sum(kx.*hk).*u(:,i-1)+sigmau*randn*sqrt(dt);
else
end
psik(1,i)=u(:,i);
%% compute vel
    xy=[x(:,i),y(:,i)];
    sx=-dY.*psik(:,i);sy=dX.*psik(:,i);
    G = exp(1i * xy * kk*2*pi/L1);
    ux=G*sx(:);uy=G*sy(:);
% % %     G1 = (exp(1i * xy * kk*2*pi/L1) .* (-ones(L,1) * transpose(dY))); % Fourier bases for u
% % %     ux1=G1*psik(:,i);rnorm(ux1(:),ux(:));

    if max(abs(imag(uy)))>10^(-6) ||  max(abs(imag(ux)))>10^(-6)
      disp('complex velocity'  )
      max(abs(imag(u)))
    else
        ux=real(ux);
        uy=real(uy);
    end

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
ind3=getind(2,0,kk);ind4=getind(2,2,kk);
ind1=getind(0,0,kk);ind2=getind(1,0,kk);
traj=real(psik([ind1,ind2,ind3],ts*gap:gap:end));
save traj2d traj

nx=30;ny=nx;Dim_Grid = nx;L1=2*pi;
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
 G1fixed = (exp(1i * xy * kk*2*pi/L1) .* (-ones(nx*ny,1) * transpose(dY))); 
 G2fixed = (exp(1i * xy * kk*2*pi/L1) .* (ones(nx*ny,1)* transpose(dX)));


return
data=real(psik(getind(2,1,kk),:));
plot(data)
plotpdf(data);

it=10000;
nx=30;ny=nx;Dim_Grid = nx;L1=2*pi;
x1d=linspace(0,L1,nx);x1d=x1d';
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
    sx=-dY.*psik(:,it);sy=dX.*psik(:,it);
    G = exp(1i * xy * kk*2*pi/L1);
    ux=G*sx(:);uy=G*sy(:);
    ux=reshape(real(ux),ny,nx);uy=reshape(real(uy),ny,nx);
    figure();
    scatter(x(:,it),y(:,it),'linewidth',1.5);hold on
quiver(xx,yy,ux,uy,1, 'linewidth',1.5,'color','r')
xlim([0,L1]);ylim([0,L1]);setgca(16)