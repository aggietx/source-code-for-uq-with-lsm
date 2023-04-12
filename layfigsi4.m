
clear;
err1=[0.4701,0.3208,0.2232,0.1610,0.1373;
    0.4757,0.3263, 0.2190,0.1621,0.1383];
err2=[ 0.4663,0.3243,0.2302,0.1773,0.1486;
    0.4853,0.3378, 0.2314,0.1764,0.1467];
figure();
subplot(1,2,1);
range=[3,6,12,24,36];
plot(range,err1(1,:),'*-k','linewidth',2);hold on;
plot(range,err1(2,:),'*-b','linewidth',2);hold on;
setgca(16);
xlabel('L','fontsize',16);
xlim([3,36]);set(gca,'XTick',range,'fontsize',16);
title('(a) Regime I','FontSize',16)

subplot(1,2,2);
range=[3,6,12,24,36];
plot(range,err2(1,:),'*-k','linewidth',2);hold on;
plot(range,err2(2,:),'*-b','linewidth',2);hold on;
setgca(16);
xlabel('L','fontsize',16);
xlim([3,36]);set(gca,'XTick',range,'fontsize',16);
title('(b) Regime II','FontSize',16)
lgnd=legend('Data-driven LSMs','Data-driven LSMs (K=7)');
set(lgnd,'FontSize',13);

% close all; %layfigsi4.fig
load  data2d_layer1_L6_p2.mat 
time1=find(uexact==max(uexact));
time2temp=find(abs(uexact-max(uexact)/2)<0.2);time2=time2temp(end-1);
time3temp=find(abs(uexact)<0.2);time3=time3temp(end);


range=1:T;
figure();
subplot(2,1,1)
plot(range,uhate(1,:),'k','linewidth',1.5);hold on
plot(range,real(estuhat(1,:)),'b','linewidth',1.5);
load  data2d_layer1_L6_p2_K7
plot(range,real(estuhat(1,:)),'r','linewidth',1.5);hold on
  setgca(16)
  
  plot(time1,real(uexact(time1)),'sk','MarkerSize',5,'linewidth',2);hold on
plot(time2,real(uexact(time2)),'sk','MarkerSize',5,'linewidth',2);hold on
plot(time3,real(uexact(time3)),'sk','MarkerSize',5,'linewidth',2);hold on

title('(a) Regime I','FontSize',16)
% xlabel('T','fontsize',14);
ylabel('u','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

xlim([200,400])
 load  data2d_layer2_L6_p2.mat 
 time1=find(uexact==max(uexact));
time2temp=find(abs(uexact-max(uexact)/2)<0.2);time2=time2temp(end-1);
time3temp=find(abs(uexact)<0.2);time3=time3temp(end);


subplot(2,1,2)
plot(range,uhate(1,:),'k','linewidth',1.5);hold on
plot(range,real(estuhat(1,:)),'b','linewidth',1.5);hold on

load  data2d_layer2_L6_p2_K7
plot(range,real(estuhat(1,:)),'r','linewidth',1.5);hold on
  plot(time1,real(uexact(time1)),'sk','MarkerSize',5,'linewidth',2);hold on
plot(time2,real(uexact(time2)),'sk','MarkerSize',5,'linewidth',2);hold on
plot(time3,real(uexact(time3)),'sk','MarkerSize',5,'linewidth',2);hold on

  setgca(16)
title('(b) Regime II','FontSize',16)
xlabel('t','fontsize',14);
ylabel('u','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
lgnd=legend('Truth','Data-driven LSMs','Data-driven LSMs (K=7)');
set(lgnd,'FontSize',13);
xlim([50,300])