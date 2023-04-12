% Ek=sum(real(psik).^2+imag(psik).^2,2);Ek=cumsum(Ek);Ek=Ek/Ek(end);
% % u5=u(1:2500/dt);u5=u(1:100:end);
load data5;load data75;load data1;
load acfv5;load acfv75;load acfv1;
load Ek1;load Ek75;load Ek5;

  
close all;
load u1;load u5;load u75;
dt=0.1;dt0=1/1000;
figure();

subplot(3,8,[1,2,3,4]);
plot(0.1:0.1:2500,u1);
xlim([0,2500])
title('(a) Trajectory of u','fontsize',14)
ylabel('p=1','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
setgca(14)
subplot(3,8,5);
plot(data1(1,:),data1(2,:),'b','linewidth',1.5);hold on;
plot(data1(3,:),data1(4,:),'--k','linewidth',1.5);
title('(b) PDF of u','fontsize',14)
setgca(14)

subplot(3,8,6);
plot(data1(1,:),log(data1(2,:)),'b','linewidth',1.5);hold on;
plot(data1(3,:),log(data1(4,:)),'--k','linewidth',1.5);
ylim([-15,1])
title('(c) PDF in log scale','fontsize',14)
setgca(14)

subplot(3,8,7);
plot(0:dt0:60,acfv1,'b','linewidth',1.5);
xlim([0,60])
setgca(14)
title('(d) ACF of u','fontsize',14)
subplot(3,8,8);
plot(1:10,Ek1,'-ob','linewidth',1.5);
xlim([1,10])
title('(e) E[\psi_{1:s}]','fontsize',14)
setgca(14)

%%%% 

subplot(3,8,[1,2,3,4]+8);
plot(0.1:0.1:2500,u75);
xlim([0,2500])
% title('(a) Trajectory of u','fontsize',14)
ylabel('p=0.75','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
setgca(14)
subplot(3,8,5+8);
plot(data75(1,:),data75(2,:),'b','linewidth',1.5);hold on;
plot(data75(3,:),data75(4,:),'--k','linewidth',1.5);
% title('(b) PDF of u','fontsize',14)
setgca(14)

subplot(3,8,6+8);
plot(data75(1,:),log(data75(2,:)),'b','linewidth',1.5);hold on;
plot(data75(3,:),log(data75(4,:)),'--k','linewidth',1.5);
ylim([-15,1])
% title('(c) PDF in log scale','fontsize',14)
setgca(14)

subplot(3,8,7+8);
plot(0:dt0:60,acfv75,'b','linewidth',1.5);
xlim([0,60])
setgca(14)
% title('(d) ACF of u','fontsize',14)
subplot(3,8,8+8);
plot(1:10,Ek75,'-ob','linewidth',1.5);
xlim([1,10])
% title('(e) E[\psi_{1:s}]','fontsize',14)
setgca(14)


%%%
subplot(3,8,[1,2,3,4]+8*2);
plot(0.1:0.1:2500,u5);
xlim([0,2500])
% title('(a) Trajectory of u','fontsize',14)
ylabel('p=0.5','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
xlabel('t','fontsize',14);
setgca(14)

subplot(3,8,5+8*2);
plot(data5(1,:),data5(2,:),'b','linewidth',1.5);hold on;
plot(data5(3,:),data5(4,:),'--k','linewidth',1.5);
% title('(b) PDF of u','fontsize',14)
setgca(14)

subplot(3,8,6+8*2);
plot(data5(1,:),log(data5(2,:)),'b','linewidth',1.5);hold on;
plot(data5(3,:),log(data5(4,:)),'--k','linewidth',1.5);
ylim([-15,1])
% title('(c) PDF in log scale','fontsize',14)
setgca(14)

subplot(3,8,7+8*2);
plot(0:dt0:60,acfv5,'b','linewidth',1.5);
xlim([0,60])
setgca(14)
% title('(d) ACF of u','fontsize',14)
subplot(3,8,8+8*2);
plot(1:10,Ek5,'-ob','linewidth',1.5);
xlim([1,10])
% title('(e) E[\psi_{1:s}]','fontsize',14)
setgca(14)
xlabel('s','fontsize',14);