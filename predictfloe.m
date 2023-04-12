clear;
% close all;
rng(1);
load  data_gb_L60_K5c;
dt0=dt;
T1=5;N1=T1/dt0;
xpredexact=zeros(L,N1);
ypredexact=zeros(L,N1);
xpredexact(:,1)=2*pi*rand(L,1);ypredexact(:,1)=2*pi*rand(L,1);
rng(5)
u_hatold=savedtraj(1:dimuhat,end);
disp('exact prediction...')
for i=2:N1
    if mod(i,floor(N1/5)) == 0;fprintf('rest step is %d\n',N1-i);end
    
    x_loc = [xpredexact(:,i-1),ypredexact(:,i-1)];
    G1 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
    G2 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
    u_hatnew=u_hatold+(alldkfit(:,end-1)).*u_hatold*dt0+sqrt(dt0)*sigmauhat*randn(dimuhat,1);
     
    u = (G1*  u_hatnew); 
    v = (G2*  u_hatnew); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity'  )
      max(abs(imag(u)))
    end
    u=real(u);
    v=real(v);
    xpredexact(:,i) = xpredexact(:,i-1) + u * dt0 + sqrt(dt0) * sig_ex * randn(L,1); % floe equation in x
    ypredexact(:,i) = ypredexact(:,i-1) + v * dt0 + sqrt(dt0) * sig_ex * randn(L,1);
    xpredexact(:,i) = mod(xpredexact(:,i),L1);
   ypredexact(:,i) = mod(ypredexact(:,i),L1);     
    u_hatold=u_hatnew;
end

xpred=zeros(L,N1);xpred(:,1)=xpredexact(:,1);
ypred=zeros(L,N1);ypred(:,1)=ypredexact(:,1);

% u_hatold=savedtraj(1:dimuhat,end);
u_hatold=savedtraj(1+2*dimuhat:3*dimuhat,end);
disp('estimated prediction...')
for i=2:N1
    if mod(i,floor(N1/5)) == 0;fprintf('rest step is %d\n',N1-i);end
    x_loc = [xpred(:,i-1),ypred(:,i-1)];
    G1 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
    G2 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
    u_hatnew=u_hatold+(-alldkfit(:,Nit)+1i*allomegakfit(:,Nit)).*u_hatold*dt0+...
         allFfit(:,Nit)*dt0+1*sqrt(dt0)*sigmafit*randn(dimuhat,1);
     
    u = (G1*  u_hatnew); 
    v = (G2*  u_hatnew); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity'  )
      max(abs(imag(u)))
    end
    u=real(u);
    v=real(v);
    xpred(:,i) = xpred(:,i-1) + u * dt0 + sqrt(dt0) * sig_ex * randn(L,1); % floe equation in x
    ypred(:,i) = ypred(:,i-1) + v * dt0 + sqrt(dt0) * sig_ex * randn(L,1);
   xpred(:,i) = mod(xpred(:,i),L1);
   ypred(:,i) = mod(ypred(:,i),L1);   
    u_hatold=u_hatnew;
end
range=1:10;
figure();
for id=1:16
subplot(4,4,id)
% id=1;
plot(xpredexact(id,range),ypredexact(id,range),'k','linewidth',1.5);hold on;
plot(xpred(id,range),ypred(id,range),'b','linewidth',1.5);hold on;
plot(xpred(id,1),ypred(id,1),'*r','linewidth',3);
% xlim([0,2*pi]);
% ylim([0,2*pi]);
setgca(14);
end
return
for i=1:10:N1
    close all;
scatter(xpredexact(:,i),ypredexact(:,i),'r')
hold on;
scatter(xpred(:,i),ypred(:,i),'b')
    title(['t = ', num2str(dt0*i)],'FontSize',24)
setgca(18)
pause()
end
t=.1;i=t/dt0;
% plot2line(xpredexact(:,i),xpred(:,i))
% plot2line(ypredexact(:,i),ypred(:,i))
figure();
plot(1:L,xpredexact(:,i),'-*r','linewidth',1.5);hold on
plot(1:L,xpred(:,i),'-ob','linewidth',1.5);
title(['x coordinate comparison at t = ', num2str(dt0*i)],'FontSize',24)
lgnd=legend('Exact','Estimated');
set(lgnd,'FontSize',18);
setgca(18)