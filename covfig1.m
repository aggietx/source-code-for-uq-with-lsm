%%% gb/total gr=2
%  ent=[ent;mean(sigfilgb(5:end)),mean(sigfilgr(5:end)),mean(sigfil(5:end)),...
%      mean(disfilgb(5:end)),mean(disfilgr(5:end)),mean(disfil(5:end))];

close all;
clear;
Lgb=[3,6,12:12:48];lim1=[min(Lgb),max(Lgb)];
Lgr=[12,30,48:36:160];lim2=[min(Lgr),max(Lgr)];
% load xexact1;load yexact1;%%% sw
% load xexact;load yexact;%%% gb 
%load dataswrand12;%%% random selection

load entswcase1;load entswcase2;load entswcase3;
load rccswcase1;load rccswcase2;load rccswcase3;%%%% gb gr full
load dataswcase1_L30;entswcase1=[entswcase1(1,:);ent;entswcase1(2:end,:)];
rccswcase1=[rccswcase1(1,:);rcc;rccswcase1(2:end,:)];
load dataswcase2_L30;entswcase2=[entswcase2(1,:);ent;entswcase2(2:end,:)];
rccswcase2=[rccswcase2(1,:);rcc;rccswcase2(2:end,:)];
load dataswcase3_L30;entswcase3=[entswcase3(1,:);ent;entswcase3(2:end,:)];
rccswcase3=[rccswcase3(1,:);rcc;rccswcase3(2:end,:)];
% return
entswcase1=[14.92  23.14 38.06 6.08 13.22 25.64
    10.33 18.93 29.27  11.02 20.49 39.76
    8.62 17.15 25.78 14.11 24.98 48.28
    6.62 15.00 21.62 18.39 31.16 59.95
    5.69 13.72 19.41 21.20 35.26 67.64
    5.11 12.92 18.02 23.27 38.31 73.38];
entswcase2=[17.15 26.81 43.96 20.10 38.31 58.42
            12.35 22.59 34.94 28.69 53.46 82.15
            10.59 21.06 31.65 33.43 62.11 95.54
            8.64 19.24 27.89 39.30 73.02 112.33
            7.62 18.07 25.70 43.15 80.28 123.43
            7.06 17.34 24.41 46.02 85.74 131.77
             ];
entswcase3=[17.20  26.15  43.35 21.93 46.87 68.81
            12.50  22.27  34.77 30.68 63.44  94.12
            10.76  20.78  31.55 35.48 72.66 108.14 
            8.83   19.01  27.84 41.41 84.15 125.56
            7.81   17.88  25.69 45.29 91.70 137.00
            7.24  17.17  24.41 48.19 97.36 145.55];

errgb=[0.69  0.72 0.80  0.63  0.76  0.65
    0.56  0.83 0.64  0.78  0.62  0.78
    0.40  0.91  0.46  0.89 0.47  0.88
    0.29  0.96  0.31  0.95 0.32  0.95
     0.24  0.97  0.26  0.97  0.26  0.96
     0.22  0.98 0.23  0.97  0.23  0.97];
% entgb=[16.47 10.24  22.43  21.91 19.72 50.95
%     10.67 16.87   14.45 28.3 13.57 50.95
%     5.81 26.21    7.70 35.24 7.87 50.95
%     3.00 37.08    3.59 42.54 3.84 50.95
%     2.16 43.09   2.47 46.95 2.82 50.95
%     1.75 47.16 1.92 50.13 1.95 50.95];

entgb=[16.23 10.31  22.00 21.91 21.83 21.93
    10.54 16.85  14.33 28.32 14.19 28.47
    5.86 26.33    7.52 35.24  7.51 35.48
    3.00 37.09    3.57 42.54 3.61 42.86
    2.16 43.14   2.45 46.95 2.48 47.30
    1.75 47.20   1.93 50.13 1.95 50.50];

% s=entgb(:,end-1);
% plot(Lgb,s,'*g','linewidth',2);

figure();

subplot(3,4,1);
plot(Lgb,errgb(:,1),'-*b','linewidth',2);hold on
plot(Lgb,errgb(:,3),'-*r','linewidth',2);hold on
plot(Lgb,errgb(:,5),'-*g','linewidth',2);hold on
setgca(13);xlim(lim1);set(gca,'XTick',Lgb,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(a) GB','FontSize',13);
subplot(3,4,5);
plot(Lgb,entgb(:,1),'-*b','linewidth',2);hold on
plot(Lgb,entgb(:,3),'-*r','linewidth',2);hold on
plot(Lgb,entgb(:,5),'-*g','linewidth',2);hold on
setgca(13);xlim(lim1);set(gca,'XTick',Lgb,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
ylabel('Signal','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)

subplot(3,4,9);
plot(Lgb,entgb(:,2),'-*b','linewidth',2);hold on
plot(Lgb,entgb(:,4),'-*r','linewidth',2);hold on
plot(Lgb,entgb(:,6),'-*g','linewidth',2);hold on
setgca(13);xlim(lim1);set(gca,'XTick',Lgb,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
xlabel('L','FontSize',14);
ylabel('Dispersion','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,4,2);
plot(Lgr,rccswcase1(:,6),'-*b','linewidth',2);hold on
plot(Lgr,rccswcase2(:,6),'-*r','linewidth',2);hold on
plot(Lgr,rccswcase3(:,6),'-*g','linewidth',2);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
% ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(b) SW','FontSize',13);

subplot(3,4,3);
plot(Lgr,rccswcase1(:,4),'-*b','linewidth',2);hold on
plot(Lgr,rccswcase2(:,4),'-*r','linewidth',2);hold on
plot(Lgr,rccswcase3(:,4),'-*g','linewidth',2);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
% ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(c) GB of SW','FontSize',13);

subplot(3,4,4);
plot(Lgr,rccswcase1(:,5),'-*b','linewidth',2);hold on
plot(Lgr,rccswcase2(:,5),'-*r','linewidth',2);hold on
plot(Lgr,rccswcase3(:,5),'-*g','linewidth',2);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
% ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(d) GR of SW','FontSize',13);
ylim([0,1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,4,6);
plot(Lgr,entswcase1(:,3),'-*b','linewidth',2);hold on
plot(Lgr,entswcase2(:,3),'-*r','linewidth',2);hold on
plot(Lgr,entswcase3(:,3),'-*g','linewidth',2);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
% ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('SW','FontSize',13);

subplot(3,4,7);
plot(Lgr,entswcase1(:,1),'-*b','linewidth',2);hold on
plot(Lgr,entswcase2(:,1),'-*r','linewidth',2);hold on
plot(Lgr,entswcase3(:,1),'-*g','linewidth',2);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
% ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('SW','FontSize',13);

subplot(3,4,8);
plot(Lgr,entswcase1(:,2),'-*b','linewidth',2);hold on
plot(Lgr,entswcase2(:,2),'-*r','linewidth',2);hold on
plot(Lgr,entswcase3(:,2),'-*g','linewidth',2);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr,'fontsize',13)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,4,6+4);
plot(Lgr,entswcase1(:,3+3),'-*b','linewidth',2);hold on
plot(Lgr,entswcase2(:,3+3),'-*r','linewidth',2);hold on
plot(Lgr,entswcase3(:,3+3),'-*g','linewidth',2);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
% ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('SW','FontSize',13);
xlabel('L','FontSize',14);
subplot(3,4,7+4);
plot(Lgr,entswcase1(:,1+3),'-*b','linewidth',2);hold on
plot(Lgr,entswcase2(:,1+3),'-*r','linewidth',2);hold on
plot(Lgr,entswcase3(:,1+3),'-*g','linewidth',2);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr,'fontsize',13)
% lgnd=legend('Full R','Diagonal R','Constant R');
% set(lgnd,'FontSize',12);
% ylabel('RMS','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('SW','FontSize',13);
xlabel('L','FontSize',14);

subplot(3,4,8+4);
plot(Lgr,entswcase1(:,2+3),'-*b','linewidth',2);hold on
plot(Lgr,entswcase2(:,2+3),'-*r','linewidth',2);hold on
plot(Lgr,entswcase3(:,2+3),'-*g','linewidth',2);hold on
setgca(13);xlim(lim2);set(gca,'XTick',Lgr,'fontsize',13)
xlabel('L','FontSize',14);
lgnd=legend('Full R','Diagonal R','Constant R');
set(lgnd,'FontSize',12);



