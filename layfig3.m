clear 
close all;
% lim=[10,110];

% load data1d_layer2_L12
load data1d_layer2_L6_p2
time1=find(uexact==max(uexact));
time2=find(abs(uexact-max(uexact)/2)<0.2);time2=time2(1);
time3=find(abs(uexact)<0.2);time3=time3(1);

lim=[1,T];
range=1:T;
% load data1d_layer2_L12
% load data1dreduced_layer2_L12
% load data_enkbf2_L12
% load data_enkbf2_reduced_L12 

figure();
subplot(3,1,1)
plot(range,real(uhate(1,:)),'k','linewidth',1);hold on
% plot(range,estuhat(1,:),'c','linewidth',1);hold on

plot(range,real(estuhat(1,:)),'b','linewidth',1.0);hold on
% load data_enkbf2_L12
load data_enkbf2_L6_p2
plot(range,real(estuhat(1,:)),'r','linewidth',1.0);hold on
% load data_enkbf2_reduced_L12 
% plot(range,estuhat(1,:),'r','linewidth',1.0);hold on

plot(time1,real(uexact(time1)),'sk','MarkerSize',5,'linewidth',2);hold on
plot(time2,real(uexact(time2)),'sk','MarkerSize',5,'linewidth',2);hold on
plot(time3,real(uexact(time3)),'sk','MarkerSize',5,'linewidth',2);hold on
% ylabel('(I) u','FontSize',16);set(get(gca,'YLabel'),'Rotation',0)
setgca(16);
% xlabel('t','fontsize',16);setgca(16);
xlim(lim)
title('Trajectory comparison');

subplot(3,1,2)
plot(range,real(uhate(2,:)),'k','linewidth',1);hold on
% plot(range,estuhat(1,:),'c','linewidth',1);hold on
% load data1d_layer2_L12
load data1d_layer2_L6_p2
plot(range,real(estuhat(2,:)),'b','linewidth',1.0);hold on
% load data_enkbf2_L12
load data_enkbf2_L6_p2
plot(range,real(estuhat(2,:)),'r','linewidth',1.0);hold on
% load data_enkbf2_reduced_L12 
% plot(range,estuhat(1,:),'r','linewidth',1.0);hold on
% ylabel('(II) \psi_1','FontSize',16);set(get(gca,'YLabel'),'Rotation',0)
setgca(16);
% xlabel('t','fontsize',16);setgca(16);
xlim(lim)


subplot(3,1,3)
plot(range,real(uhate(4,:)),'k','linewidth',1);hold on
% plot(range,estuhat(1,:),'c','linewidth',1);hold on
% load data1d_layer2_L12
 load data1d_layer2_L6_p2
plot(range,real(estuhat(4,:)),'b','linewidth',1.0);hold on
% load data_enkbf2_L12
load data_enkbf2_L6_p2
plot(range,real(estuhat(4,:)),'r','linewidth',1.0);hold on
% load data_enkbf2_reduced_L12 
% plot(range,estuhat(1,:),'r','linewidth',1.0);hold on
% ylabel('(III) \psi_2','FontSize',16);set(get(gca,'YLabel'),'Rotation',0)

xlabel('t','fontsize',16);setgca(16);
xlim(lim)
lgnd=legend('Truth','1D estimation ','EnKBF all modes');
set(lgnd,'FontSize',16);





