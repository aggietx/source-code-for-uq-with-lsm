clear;
% close all
load trajc2;load ufreec2
% load trajc1;load ufreec1
dt=0.01;ts=10;Naut=1000;
range=ts/dt:5000/dt;
figure()
for i=1:7
    subplot(7,2,2*i-1)
    [prob,r1]=ksdensity(real(traj(i,range)));
    plot(r1,prob,'k','Linewidth',2);
    [prob,r1]=ksdensity(real(ufree(i,range)));hold on
    plot(r1,prob,'b','Linewidth',2);
    setgca(16)
    
    subplot(7,2,2*i)
    if i==1
    acf1=autocorr(real(traj(i,range)),Naut*10);
    plot(0:dt:Naut*10*dt,acf1,'k','Linewidth',2);
    acf1=autocorr(real(ufree(i,range)),Naut*10);hold on
    plot(0:dt:Naut*10*dt,acf1,'b','Linewidth',2);    
    else
     acf1=autocorr(real(traj(i,range)),Naut);
    plot(0:dt:Naut*dt,acf1,'k','Linewidth',2); 
    acf1=autocorr(real(ufree(i,range)),Naut);hold on
    plot(0:dt:Naut*dt,acf1,'b','Linewidth',2);    
    end    
    setgca(16)
end
  lgnd=legend('Truth',...
    'Free run fit');
set(lgnd,'FontSize',12);
return
traj=traj(:,10:10:end);save trajc2 traj
ufree=ufree(:,10:10:end);save ufreec2 ufree
%  sampleexact=psik(:,10000:10000:50000);save  sampleexact  sampleexact
%  sampleest=gamma_mean_trace(:,10000:10000:50000);save  sampleest  sampleest
% get sampleest.mat sampleexact.mat
load sampleexact;load sampleest
 ii=4;
 ux=reshape(real(G1fixed* sampleexact(:,ii)),ny,nx)';
  uy=reshape(real(G2fixed* sampleexact(:,ii)),ny,nx)';

 uestx=reshape(real(G1fixed* sampleest(:,ii)),ny,nx)';
  uesty=reshape(real(G2fixed* sampleest(:,ii)),ny,nx)';
  figure();
subplot(1,2,1)
quiver(xx,yy,ux,uy);xlim([0,pi*2]);ylim([0,pi*2]);setgca(18)
subplot(1,2,2)
quiver(xx,yy,uestx,uesty)
xlim([0,pi*2]);ylim([0,pi*2]);setgca(18)
return
 ux1=reshape(real(G1fixed* psik(:,40000)),ny,nx)';
  uy1=reshape(real(G2fixed*psik(:,40000)),ny,nx)';
figure();quiver(xx,yy,ux1,uy1)
xlim([0,pi*2]);ylim([0,pi*2]);setgca(18)
return
clear
load  data_gb_L12_K5c;
s=savedtraj(1:dimuhat,10);
u=ifft2(reshape(s,11,11));
 xy=kk'*2*pi/(2*K_max+1);
G1 = (exp(1i * xy * kk))/(2*K_max+1)^2;
u1=G1*s(:);
u1=reshape(u1,11,11);
return
computeGBvel(rk,kk,L1,Dim_Grid,uhat);

return
ns=2*pi;Dim_Grid=2*K_max+1;
[xx,yy] = meshgrid(linspace(0,ns,Dim_Grid), linspace(0,ns,Dim_Grid));
    x_loc = [xx(:),yy(:)];

 G =(exp(1i * x_loc * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * ones(1,Dim_Grid^2)));
 iG =(exp(-1i * x_loc * kk*2*pi/L1).* (ones(Dim_Grid^2,1) * ones(1,Dim_Grid^2)));
 
return
ind1=getind(2,2,kk);inda=ind1(1);
rmscc(real(gamma_mean_trace(inda,:)),real(u_hat(inda,:)),1);
rmscc(real(gamma_mean_trace_smoother(inda,:)),real(u_hat(inda,:)),1);

return
a=real(u_hat);
b=imag(u_hat);
c=sum(a.^2,2)+sum(b.^2,2);
sum(c(dimuhat0*1+1:end))/sum(c(1:dimuhat0))
Dim_Grid = 50;
ns=L1;
[xx,yy] = meshgrid(linspace(0,ns,Dim_Grid), linspace(0,ns,Dim_Grid));
x_vec = [reshape(xx,[],1), reshape(yy,[],1)]; 
    x_loc = [xx(:),yy(:)];
    G1 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(1,:))); % Fourier bases for u
    G2 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(2,:))); % Fourier bases for v

for i=20:10:N
% for i=15000
    u = (G1*  u_hat(:,i-1)); 
    v = (G2*  u_hat(:,i-1)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity'  )
    end
    u=real(u);
    v=real(v);
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);

    imagesc([0,ns],[0,ns],sqrt(real(u).^2+real(v).^2));colorbar;setgca(16);
    hold on
    quiver(xx, yy, u, v, 'linewidth',1)
% hold on;
% 
    xlim([0, ns ])
    ylim([0, ns ])    
%     xlim([-pi, pi ])
%     ylim([-pi, pi ])
    box on    
    title(['t = ', num2str(i*dt)],'fontsize',16)
%     pause(.001);
hold on
scatter(xexact(1,i),yexact(1,i),'r')
pause()
% close all
%     u0=[u0,sqrt(u(:).^2+v(:).^2)];
    setgca(16);
end
