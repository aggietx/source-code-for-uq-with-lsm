clear;
% close all;
load rmsregII1;load rmsregII2;load rmsregII3;load rmsregII4;load rmsregII5;
load rmsregI1;load rmsregI2;load rmsregI3;load rmsregI4;load rmsregI5;
LL=[3,6,12,24,36];
xli=[3,36];
figure()
for i=1:3;

subplot(2,3,i);
plot(LL,rmsregI1(:,i),'-*k','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
if i==1;ylabel('(I) Regime I','FontSize',16);set(get(gca,'YLabel'),'Rotation',0);end

plot(LL,rmsregI2(:,i),'-*g','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
plot(LL,rmsregI3(:,i),'-*b','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
plot(LL,rmsregI4(:,i),'-*c','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
plot(LL,rmsregI5(:,i),'-*r','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
if i==1
    title('RMS of u','fontsize',16);
elseif i==2
    title('RMS of 2 modes','fontsize',16);
else
    title('RMS of all modes','fontsize',16);
end


subplot(2,3,i+3);
plot(LL,rmsregII1(:,i),'-*k','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
if i==1;ylabel('(I) Regime II','FontSize',16);set(get(gca,'YLabel'),'Rotation',0);end
plot(LL,rmsregII2(:,i),'-*g','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
plot(LL,rmsregII3(:,i),'-*b','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
plot(LL,rmsregII4(:,i),'-*c','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
plot(LL,rmsregII5(:,i),'-*r','linewidth',1.5);hold on;setgca(14);xlim(xli);set(gca,'XTick',LL,'fontsize',16)
if i==3
lgnd=legend('EnKBF with all modes','EnKBF with 2 modes','1D estimation with all modes','1D estimation with 2 modes','2D estimation');
set(lgnd,'FontSize',14);
end
xlabel('L','fontsize',14);
end