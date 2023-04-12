clear;
close all;
m=3;%%%% ensemble size
gamma=0.25;
T=5000;
dt=0.01;
N=T/dt;
sigma=sqrt(0.35);
q=zeros(N,1);
v=zeros(N,1);
Q=zeros(N,1);
c=0.5;
c1=0.5;%%% initial error
q(1)=1;
v(1)=1;
rng(1);
for i=1:N-1
    q(i+1)=q(i)+v(i)*dt;
    v(i+1)=v(i)-(-sin(q(i))+0.5*(q(i)/6)^3+0.1 )*dt-gamma*v(i)*dt+sigma*randn*sqrt(dt);
    Q(i+1)=Q(i)+v(i)*dt+sqrt(c)*sqrt(dt)*randn;
end
dQ=Q(2:end)-Q(1:end-1);
figure();plot(q);
% figure();plot(Q);

vens=zeros(m,N);
qens=zeros(m,N);
vens(:,1)=v(1)*ones(m,1)+c1*randn(m,1);
qens(:,1)=q(1)*ones(m,1)+c1*randn(m,1);
vens(:,2)=v(2)*ones(m,1)+c1*randn(m,1);
qens(:,2)=q(2)*ones(m,1)+c1*randn(m,1);

for i=2:N-1
    vbar=mean(vens(:,i));qbar=mean(qens(:,i));
    P_qv=(qens(:,i)-qbar)'*(vens(:,i)-vbar)/(m-1);
    P_vv=(vens(:,i)-vbar)'*(vens(:,i)-vbar)/(m-1);
    qens(:,i+1)=qens(:,i)+vens(:,i)*dt-P_qv/2/c*( vens(:,i)*dt+vbar*dt-2*(Q(i)-Q(i-1)));
    vens(:,i+1)=vens(:,i)-(-sin(qens(:,i))+0.5*(qens(:,i)/6).^3+0.1 )*dt-...
        gamma*vens(:,i)*dt+sigma*randn(m,1)*sqrt(dt)-P_vv/2/c*( vens(:,i)*dt+vbar*dt-2*(Q(i)-Q(i-1)));
end
% figure();plot(mean(qens));
% plot2line(q,mean(qens));
subplot(2,1,1)
plot(dt:dt:T,q,'r',...
    dt:dt:T,mean(qens),'g');
title(['c is ', num2str(c), ',  ens is ', num2str(m)],'FontSize',16) 
lgnd=legend('true'  ,  'posterior');
set(lgnd,'FontSize',16,'Location', 'Best');
set(gca,'FontSize',16)

subplot(2,1,2);plot(dt:dt:T,q-mean(qens)','r',...
    dt:dt:T,mean(qens),'g');
title(['c is ', num2str(c), ',  ens is ', num2str(m)],'FontSize',16) 
lgnd=legend('diff'  ,  'posterior');
set(lgnd,'FontSize',16,'Location', 'Best');
set(gca,'FontSize',16)

return
%% pure forcast
q1=zeros(N,1);
v1=zeros(N,1);
q1(1)=1+c1*randn;
v1(1)=1+c1*randn;
rng(1);
for i=1:N-1
    q1(i+1)=q1(i)+v1(i)*dt;
    v1(i+1)=v1(i)-(-sin(q1(i))+0.5*(q1(i)/6)^3+0.1 )*dt-gamma*v1(i)*dt+sigma*randn*sqrt(dt);
end
figure();
plot(dt:dt:T,q,'r',...
    dt:dt:T,q1,'g');
lgnd=legend('true'  ,  'forcast');
set(lgnd,'FontSize',16,'Location', 'Best');
set(gca,'FontSize',16)
return
%% enkf

vens=zeros(m,N);
qens=zeros(m,N);
Qens=zeros(m,N);
vens(:,1)=v(1)*ones(m,1)+c1*randn(m,1);
qens(:,1)=q(1)*ones(m,1)+c1*randn(m,1);
Qens(:,1)
for i=1:N-1
    %% forcast
    qens(:,i+1)=qens(:,i)+vens(:,i)*dt;
    vens(:,i+1)=vens(:,i)-(-sin(qens(:,i))+0.5*(qens(:,i)/6).^3+0.1 )*dt-...
        gamma*vens(:,i)*dt+sigma*randn(m,1)*sqrt(dt);
    Qens(:,i+1)=Qens(:,i)+vens(:,i)*dt+sqrt(dt)*randn(m,1)*sqrt(dt);
    %% analysis
    vbar=mean(vens(:,i));qbar=mean(qens(:,i));
    vbar=mean(vens(:,i+1));qbar=mean(qens(:,i+1));
    P_qv=(qens(:,i+1)-qbar)'*(vens(:,i+1)-vbar)/(m-1);
    P_vv=(vens(:,i+1)-vbar)'*(vens(:,i+1)-vbar)/(m-1);    
end
