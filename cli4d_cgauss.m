%%% generate 4d climiate model
clear;
close all;
rng(3);
L12=1;L13=0.5;L24=0.5;a1=2;a2=1;d1=-1;d2=-.4;eps=1;
sigma1=0.5;sigma2=2;sigma3=0.5;sigma4=1;b123=1.5;b213=1.5;b312=-b213-b213;F1=0;F2=0;F3=0;F4=0;
T=100;dt=1/50;N=T/dt;Dt=0.1;sigma_obs=0.1;dis_obs=Dt/dt;
x1=zeros(1,N);kf=1;%%% kalman filter or not
x2=x1;x3=x1;x4=x1;
gamma1=1;gamma2=1;
MM=1;MM1=50;r=0.1;
for i=2:N
    x1(:,i)=x1(:,i-1)+dt*( -x2(:,i-1)*(L12+a1*x1(:,i-1)+a2*x2(:,i-1))+d1*x1(:,i-1)+F1*ones(MM,1)+...
    +L13*x3(:,i-1)+b123*x2(:,i-1)*x3(:,i-1))+sigma1*sqrt(dt)*randn(MM,1);
    x2(:,i)=x2(:,i-1)+dt*( x1(:,i-1)*(L12+a1*x1(:,i-1)+a2*x2(:,i-1))+d2*x2(:,i-1)+F2*ones(MM,1)+...
    +L24*x4(:,i-1)+b213*x1(:,i-1)*x3(:,i-1))+sigma2*sqrt(dt)*randn(MM,1);  
    x3(:,i)=x3(:,i-1)+dt*(-L13*x1(:,i-1)+b312*x1(:,i-1)*x2(:,i-1)+F3*ones(MM,1)-...
        gamma1/eps*x3(:,i-1))+sigma3/sqrt(eps)*sqrt(dt)*randn(MM,1);  
    x4(:,i)=x4(:,i-1)+dt*(-L24*x2(:,i-1)+F4*ones(MM,1)-gamma2/eps*x4(:,i-1))+...
        sigma4/sqrt(eps)*sqrt(dt)*randn(MM,1);      
end
% plotpdf(x1);
figure();
subplot(4,3,[1,2])
plot(dt:dt:T,x1,'linewidth',1.5);set(gca,'FontSize',15);
subplot(4,3,[4,5])
plot(dt:dt:T,x2,'linewidth',1.5);set(gca,'FontSize',15);
subplot(4,3,[7,8])
plot(dt:dt:T,x3,'linewidth',1.5);set(gca,'FontSize',15);
subplot(4,3,[10,11])
plot(dt:dt:T,x4,'linewidth',1.5);set(gca,'FontSize',15);

subplot(4,3,3);
[prob,xp]=ksdensity(x1);
plot(xp,prob,'linewidth',1.5);set(gca,'FontSize',15);
subplot(4,3,6);
[prob,xp]=ksdensity(x2);
plot(xp,prob,'linewidth',1.5);set(gca,'FontSize',15);
subplot(4,3,9);
[prob,xp]=ksdensity(x3);
plot(xp,prob,'linewidth',1.5);set(gca,'FontSize',15);
subplot(4,3,12);
[prob,xp]=ksdensity(x4);
plot(xp,prob,'linewidth',1.5);set(gca,'FontSize',15);

%% kalman filter
MM=MM1;
x1old=x1(:,1)*ones(MM,1)+sigma_obs*ones(MM,1);
x2old=x2(:,1)*ones(MM,1)+sigma_obs*ones(MM,1);
x3old=x3(:,1)*ones(MM,1)+sigma_obs*ones(MM,1);
x4old=x4(:,1)*ones(MM,1)+sigma_obs*ones(MM,1);
allpri=zeros(N/Dt,4,MM);
allpost=zeros(N/Dt,4,MM);
 M_cov=sigma_obs.^2*eye(2);
 Ro=sigma_obs.^2*speye(2);
 obsx1=x1+sigma_obs*ones(1,N);
 obsx2=x2+sigma_obs*ones(1,N);
 xpostmean=zeros(4,N/dis_obs);
H=speye(2,4);


for i=2:N
    
    
end