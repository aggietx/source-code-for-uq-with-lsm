clear;
close all

Lgb=[1:6,12:6:48];lim1=[min(Lgb),max(Lgb)];
Lgr=[1:2:12,24:24:144];lim2=[min(Lgr),max(Lgr)];
subplot(3,4,1);
load datagbrandLL
plot(Lgb,rcc(:,2),'-*g','linewidth',1);hold on
load datagbrand12
plot(Lgb(1:length(rcc)),rcc(:,2),'-*c','linewidth',1);hold on
load datagbrand24
plot(Lgb(1:length(rcc)),rcc(:,2),'-*r','linewidth',1);hold on
load datagbrand36
plot(Lgb(1:length(rcc)),rcc(:,2),'-*b','linewidth',1);hold on
load datagbrand48
plot(Lgb(1:length(rcc)),rcc(:,2),'-*k','linewidth',1);hold on

setgca(13);xlim(lim1);set(gca,'XTick',Lgb([1,6,7:end]),'fontsize',13)
ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(a) GB','FontSize',13);

subplot(3,4,5)
load datagbrandLL;ent(7,1)=5.37;
plot(Lgb,ent(:,1),'-*g','linewidth',1);hold on
load datagbrand12
plot(Lgb(1:length(ent)),ent(:,1),'-*c','linewidth',1);hold on
load datagbrand24
plot(Lgb(1:length(ent)),ent(:,1),'-*r','linewidth',1);hold on
load datagbrand36
plot(Lgb(1:length(ent)),ent(:,1),'-*b','linewidth',1);hold on
load datagbrand48
plot(Lgb(1:length(ent)),ent(:,1),'-*k','linewidth',1);hold on

setgca(13);xlim(lim1);set(gca,'XTick',Lgb([1,6,7:end]),'fontsize',13)
ylabel('Signal','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('(a) GB','FontSize',13);

subplot(3,4,9)
load datagbrandLL;
plot(Lgb,ent(:,2),'-*g','linewidth',1);hold on
load datagbrand12
plot(Lgb(1:length(ent)),ent(:,2),'-*c','linewidth',1);hold on
load datagbrand24
plot(Lgb(1:length(ent)),ent(:,2),'-*r','linewidth',1);hold on
load datagbrand36
plot(Lgb(1:length(ent)),ent(:,2),'-*b','linewidth',1);hold on
load datagbrand48
plot(Lgb(1:length(ent)),ent(:,2),'-*k','linewidth',1);hold on
xlabel('L','FontSize',14);
lgnd=legend('L/L','L/12','L/24','L/36','L/48');
set(lgnd,'FontSize',12);


setgca(13);xlim(lim1);set(gca,'XTick',Lgb([1,6,7:end]),'fontsize',13)
ylabel('Dispersion','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('(a) GB','FontSize',13);

%%%%%%%%%%%%%%%%%%%%%
subplot(3,4,2)
load dataswrandLL;%%% gb gr fil
plot(Lgr(1:length(ent)),rcc(:,6),'-*g','linewidth',1);hold on
% load dataswrand11;
load dataswrand11b;
plot(Lgr(1:length(ent)),rcc(:,6),'-*c','linewidth',1);hold on
load dataswrand48;
plot(Lgr(1:length(ent)),rcc(:,6),'-*r','linewidth',1);hold on
load dataswrand96;
plot(Lgr(1:length(ent)),rcc(:,6),'-*b','linewidth',1);hold on
% load dataswrand120;
load dataswrand144;
plot(Lgr(1:length(ent)),rcc(:,6),'-*k','linewidth',1);hold on

setgca(13);xlim(lim2);set(gca,'XTick',Lgr([3,7:end]),'fontsize',13)
title('(b) SW','FontSize',13);

subplot(3,4,3)
load dataswrandLL;%%% gb gr fil
plot(Lgr(1:length(ent)),rcc(:,4),'-*g','linewidth',1);hold on
% load dataswrand11;
load dataswrand11b;
plot(Lgr(1:length(ent)),rcc(:,4),'-*c','linewidth',1);hold on
load dataswrand48;
plot(Lgr(1:length(ent)),rcc(:,4),'-*r','linewidth',1);hold on
load dataswrand96;
plot(Lgr(1:length(ent)),rcc(:,4),'-*b','linewidth',1);hold on
% load dataswrand120;
load dataswrand144;
plot(Lgr(1:length(ent)),rcc(:,4),'-*k','linewidth',1);hold on

setgca(13);xlim(lim2);set(gca,'XTick',Lgr([3,7:end]),'fontsize',13)
title('(c) GB of SW','FontSize',13);

subplot(3,4,4)
load dataswrandLL;%%% gb gr fil
plot(Lgr(1:length(ent)),rcc(:,5),'-*g','linewidth',1);hold on
% load dataswrand11;
load dataswrand11b;
plot(Lgr(1:length(ent)),rcc(:,5),'-*c','linewidth',1);hold on
load dataswrand48;
plot(Lgr(1:length(ent)),rcc(:,5),'-*r','linewidth',1);hold on
load dataswrand96;
plot(Lgr(1:length(ent)),rcc(:,5),'-*b','linewidth',1);hold on
% load dataswrand120;
load dataswrand144;
plot(Lgr(1:length(ent)),rcc(:,5),'-*k','linewidth',1);hold on

setgca(13);xlim(lim2);set(gca,'XTick',Lgr([3,7:end]),'fontsize',13)
title('(d) GB of SW','FontSize',13);
xlabel('L','FontSize',14);
%%%%%%

subplot(3,4,6)
load dataswrandLL;%%% gb gr fil
plot(Lgr(1:length(ent)),ent(:,3),'-*g','linewidth',1);hold on
% load dataswrand11;
load dataswrand11b;
plot(Lgr(1:length(ent)),ent(:,3),'-*c','linewidth',1);hold on
load dataswrand48;
plot(Lgr(1:length(ent)),ent(:,3),'-*r','linewidth',1);hold on
load dataswrand96;
plot(Lgr(1:length(ent)),ent(:,3),'-*b','linewidth',1);hold on
% load dataswrand120;
load dataswrand144;
plot(Lgr(1:length(ent)),ent(:,3),'-*k','linewidth',1);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr([3,7:end]),'fontsize',13)

subplot(3,4,7)
load dataswrandLL;%%% gb gr fil
plot(Lgr(1:length(ent)),ent(:,1),'-*g','linewidth',1);hold on
% load dataswrand11;
load dataswrand11b;
plot(Lgr(1:length(ent)),ent(:,1),'-*c','linewidth',1);hold on
load dataswrand48;
plot(Lgr(1:length(ent)),ent(:,1),'-*r','linewidth',1);hold on
load dataswrand96;
plot(Lgr(1:length(ent)),ent(:,1),'-*b','linewidth',1);hold on
% load dataswrand120;
load dataswrand144;
plot(Lgr(1:length(ent)),ent(:,1),'-*k','linewidth',1);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr([3,7:end]),'fontsize',13)

subplot(3,4,8)
load dataswrandLL;%%% gb gr fil
plot(Lgr(1:length(ent)),ent(:,2),'-*g','linewidth',1);hold on
% load dataswrand11;
load dataswrand11b;
plot(Lgr(1:length(ent)),ent(:,2),'-*c','linewidth',1);hold on
load dataswrand48;
plot(Lgr(1:length(ent)),ent(:,2),'-*r','linewidth',1);hold on
load dataswrand96;
plot(Lgr(1:length(ent)),ent(:,2),'-*b','linewidth',1);hold on
% load dataswrand120;
load dataswrand144;
plot(Lgr(1:length(ent)),ent(:,2),'-*k','linewidth',1);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr([3,7:end]),'fontsize',13)
% xlabel('L','FontSize',14);
%%%%%%
subplot(3,4,10)
load dataswrandLL;%%% gb gr fil
plot(Lgr(1:length(ent)),ent(:,3+3),'-*g','linewidth',1);hold on
% load dataswrand11;
load dataswrand11b;
plot(Lgr(1:length(ent)),ent(:,3+3),'-*c','linewidth',1);hold on
load dataswrand48;
plot(Lgr(1:length(ent)),ent(:,3+3),'-*r','linewidth',1);hold on
load dataswrand96;
plot(Lgr(1:length(ent)),ent(:,3+3),'-*b','linewidth',1);hold on
% load dataswrand120;
load dataswrand144;
plot(Lgr(1:length(ent)),ent(:,3+3),'-*k','linewidth',1);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr([3,7:end]),'fontsize',13)
xlabel('L','FontSize',14);

subplot(3,4,11)
load dataswrandLL;%%% gb gr fil
plot(Lgr(1:length(ent)),ent(:,1+3),'-*g','linewidth',1);hold on
% load dataswrand11;
load dataswrand11b;
plot(Lgr(1:length(ent)),ent(:,1+3),'-*c','linewidth',1);hold on
load dataswrand48;
plot(Lgr(1:length(ent)),ent(:,1+3),'-*r','linewidth',1);hold on
load dataswrand96;
plot(Lgr(1:length(ent)),ent(:,1+3),'-*b','linewidth',1);hold on
% load dataswrand120;
load dataswrand144;
plot(Lgr(1:length(ent)),ent(:,1+3),'-*k','linewidth',1);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr([3,7:end]),'fontsize',13)
xlabel('L','FontSize',14);

subplot(3,4,12)
load dataswrandLL;%%% gb gr fil
plot(Lgr(1:length(ent)),ent(:,2+3),'-*g','linewidth',1);hold on
% load dataswrand11;
load dataswrand11b;
plot(Lgr(1:length(ent)),ent(:,2+3),'-*c','linewidth',1);hold on
load dataswrand48;
plot(Lgr(1:length(ent)),ent(:,2+3),'-*r','linewidth',1);hold on
load dataswrand96;
plot(Lgr(1:length(ent)),ent(:,2+3),'-*b','linewidth',1);hold on
% load dataswrand120;
load dataswrand144;
plot(Lgr(1:length(ent)),ent(:,2+3),'-*k','linewidth',1);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr([3,7:end]),'fontsize',13)
xlabel('L','FontSize',14);
lgnd=legend('L/L','L/11','L/48','L/96','L/144');
set(lgnd,'FontSize',12);
