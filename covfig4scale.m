clear;
close all
%load Rdiag48;s1=mean(Rdiag48(2:end,:),1);corr2(1./sqrt(Lgb),s1);
Lgb=[1:6,12:6:48];lim1=[min(Lgb),max(Lgb)];
Lgr=[1:2:12,24:24:144];lim2=[min(Lgr),max(Lgr)];
subplot(3,1,1);
load datagbrandLLscale
plot(Lgb,rcc(:,2),'-*g','linewidth',1);hold on
load datagbrand12scale
plot(Lgb(1:length(rcc)),rcc(:,2),'-*c','linewidth',1);hold on
load datagbrand24scale
plot(Lgb(1:length(rcc)),rcc(:,2),'-*r','linewidth',1);hold on
load datagbrand36scale
plot(Lgb(1:length(rcc)),rcc(:,2),'-*b','linewidth',1);hold on
load datagbrand48scale
plot(Lgb(1:length(rcc)),rcc(:,2),'-*k','linewidth',1);hold on

setgca(13);xlim(lim1);set(gca,'XTick',Lgb([1,6,7:end]),'fontsize',13)
ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(a) GB','FontSize',13);

subplot(3,1,2)
load datagbrandLLscale;ent(7,1)=5.37;
plot(Lgb,ent(:,1),'-*g','linewidth',1);hold on
load datagbrand12scale
plot(Lgb(1:length(ent)),ent(:,1),'-*c','linewidth',1);hold on
load datagbrand24scale
plot(Lgb(1:length(ent)),ent(:,1),'-*r','linewidth',1);hold on
load datagbrand36scale
plot(Lgb(1:length(ent)),ent(:,1),'-*b','linewidth',1);hold on
load datagbrand48scale
plot(Lgb(1:length(ent)),ent(:,1),'-*k','linewidth',1);hold on

setgca(13);xlim(lim1);set(gca,'XTick',Lgb([1,6,7:end]),'fontsize',13)
ylabel('Signal','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('(a) GB','FontSize',13);

subplot(3,1,3)
load datagbrandLLscale;
plot(Lgb,ent(:,2),'-*g','linewidth',1);hold on
load datagbrand12scale
plot(Lgb(1:length(ent)),ent(:,2),'-*c','linewidth',1);hold on
load datagbrand24scale
plot(Lgb(1:length(ent)),ent(:,2),'-*r','linewidth',1);hold on
load datagbrand36scale
plot(Lgb(1:length(ent)),ent(:,2),'-*b','linewidth',1);hold on
load datagbrand48scale
plot(Lgb(1:length(ent)),ent(:,2),'-*k','linewidth',1);hold on
xlabel('L','FontSize',14);
lgnd=legend('L/L','L/12','L/24','L/36','L/48');
set(lgnd,'FontSize',12);


setgca(13);xlim(lim1);set(gca,'XTick',Lgb([1,6,7:end]),'fontsize',13)
ylabel('Dispersion','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('(a) GB','FontSize',13);