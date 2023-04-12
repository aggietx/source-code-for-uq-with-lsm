%%%% record gb data
close all;
xrange=[6,12,24,48,72,96];
Kmax=5; 
L=48;
% edk=[0.929,0.692,0.526,0.497,0.524,0.557;
%     0.929,0.582,0.376,0.334,0.349,0.374;
%     0.929,0.497,0.287,0.243,0.236,0.235;
%     0.929,0.385,0.233,0.225,0.224,0.223;
%      0.929,0.343,0.222,0.225,0.227,0.227;
%       0.929,0.309,0.208,0.212,0.212,0.212;
%  ];
% esigma=[0.843,0.448,0.365,0.368,0.383,0.398;
%     0.843,0.394,0.297,0.288,0.292,0.298;
%     0.843,0.340,0.232,0.215,0.212,0.212;
%     0.843,0.276,0.179,0.168,0.166,0.166;
%        0.843,0.266,0.192,0.184,0.184,0.184;
%           0.843,0.225,0.179,0.175,0.174,0.174];
edk=[0.9291    0.6815    0.5911    0.6280    0.6792    0.7176;
     0.9291    0.5946    0.4652    0.4657    0.5061    0.5367;
     0.9291    0.4902    0.4246    0.4399    0.4425    0.4417;
     0.9291    0.3933    0.3737    0.4103    0.4192    0.4214;
     0.9291    0.3763    0.2748    0.2760    0.2768    0.2769;
     0.9291    0.3614    0.2868    0.2825    0.2810    0.2805];
esigma=[0.8429    0.5026    0.4155    0.4225    0.4372    0.4501;
    0.8429    0.4306    0.3196    0.3117    0.3184    0.3247;
    0.8429    0.3511    0.2233    0.2103    0.2079    0.2070;
    0.8429    0.2874    0.1697    0.1634    0.1629    0.1627;
    0.8429    0.2632    0.1608    0.1527    0.1515    0.1512;
    0.8429    0.2336    0.1382    0.1323    0.1313    0.1311];
    
      figure()
subplot(2,1,1)
plot(xrange,edk(:,1),'*-b','linewidth',2);hold on;
plot(xrange,edk(:,2),'*-r','linewidth',2);hold on;
plot(xrange,edk(:,4),'*-g','linewidth',2);hold on;
plot(xrange,edk(:,6),'*-k','linewidth',2);hold on;
% lgnd=legend('Initial','1 iteration','3 iterations','5 iterations');
% set(lgnd,'FontSize',12);
setgca(18)
title('Relative L_2 error of d_k','FontSize',18);
xlim([6,96])
set(gca,'XTick',[6,12,24,48,72,96])
% xticklabels({'1','20','100','100','1000'})

subplot(2,1,2)
plot(xrange,esigma(:,1),'*-b','linewidth',2);hold on;
plot(xrange,esigma(:,2),'*-r','linewidth',2);hold on;
plot(xrange,esigma(:,4),'*-g','linewidth',2);hold on;
plot(xrange,esigma(:,6),'*-k','linewidth',2);hold on;
lgnd=legend('Initial guess','1 iteration','3 iterations','5 iterations');
set(lgnd,'FontSize',14);
setgca(18)
title('Relative L_2 error of \sigma_k','FontSize',18);
xlabel('L','FontSize',18);

xlim([6,96])
set(gca,'XTick',[6,12,24,48,72,96])
% xticklabels({'1','20','100','100','1000'})

%% spectral
load allene6;allene6(:,1:end-1)=-allene6(:,1:end-1);
load allene12;allene12(:,1:end-1)=-allene12(:,1:end-1);
load allene24;allene24(:,1:end-1)=-allene24(:,1:end-1);
load allene48;allene48(:,1:end-1)=-allene48(:,1:end-1);
load allene72;allene72(:,1:end-1)=-allene72(:,1:end-1);
load allene96;allene96(:,1:end-1)=-allene96(:,1:end-1);
      figure()
plot(allene12(:,end),'*-b','linewidth',2);hold on;
plot(allene6(:,1+Nit),'*-r','linewidth',2);hold on;
% plot(allene12(:,1+Nit),'*-k','linewidth',2);hold on;
plot(allene24(:,1+Nit),'*-g','linewidth',2);hold on;
plot(allene48(:,1+Nit),'*-k','linewidth',2);hold on;
% plot(allene72(:,1+Nit),'*-g','linewidth',2);hold on;
plot(allene96(:,1+Nit),'*-c','linewidth',2);hold on;
% plot(edk(:,6),'*-k','linewidth',2);hold on;
lgnd=legend('Exact energy','L=6','L=24','L=48','L=96');
set(lgnd,'FontSize',16);
setgca(18)
xlabel('|k|','FontSize',18);
xlim([1,size(allene12,1)]);
title('energy','FontSize',18);
% xlim([6,96])
% set(gca,'XTick',[6,12,24,48,72,96])

%% relative enetropy
extf=[120.653 382.345;
    120.782 672.909;
    121.630 1128.295;
    121.787 1782.680;
    122.007 2553.548;
    121.860 2996.514];
    
estf=[ 339.136 657.349;
    220.503 839.522;
    154.737 1204.425;
     135.903 1817.760;
     130.933 2422.522;
     127.621 2842.848];

exft=[63.366 64.644;
    44.032 98.627;
    27.794 141.886;
    16.412 189.380;
    12.148 220.635;
    9.814 240.133];
esft=[63.366 118.904;
    44.032 131.746;
    27.794 156.754;
    16.412 197.031;
    12.148 224.465;
    9.814 242.124];

      figure()
subplot(2,1,1)
plot(xrange,extf(:,1),'*-b','linewidth',2);hold on;
plot(xrange,estf(:,1),'*-r','linewidth',2);hold on;

title('Signal (true vs filter)','FontSize',18);
setgca(14)
xlabel('L','FontSize',18);
set(gca,'XTick',[6,12,24,48,72,96])
xlim([6,96])

subplot(2,1,2)
plot(xrange,extf(:,2),'*-b','linewidth',2);hold on;
plot(xrange,estf(:,2),'*-r','linewidth',2);hold on;



title('Dispersion (true vs filter)','FontSize',18);
setgca(18)
xlabel('L','FontSize',18);
set(gca,'XTick',[6,12,24,48,72,96])
xlim([6,96])
lgnd=legend('True','Estimation');
set(lgnd,'FontSize',16);

figure();
subplot(2,1,1)
plot(xrange,exft(:,1),'*-b','linewidth',2);hold on;
plot(xrange,esft(:,1),'*-r','linewidth',2);hold on;

title('Signal (filter vs true)','FontSize',18);
setgca(14)
xlabel('L','FontSize',18);
set(gca,'XTick',[6,12,24,48,72,96])
xlim([6,96])

subplot(2,1,2)
plot(xrange,exft(:,2),'*-b','linewidth',2);hold on;
plot(xrange,esft(:,2),'*-r','linewidth',2);hold on;



title('dispersion (filter vs true)','FontSize',18);
setgca(18)
xlabel('L','FontSize',18);
set(gca,'XTick',[6,12,24,48,72,96])
xlim([6,96])


lgnd=legend('True','Estimation');
set(lgnd,'FontSize',16);


%% plot estimated parameter comparison with using full R and constant R
load alldkfit;load allsigmak2;load allomegakfit %%% must reload the data
figure()
subplot(2,1,1);
plot(-alldkfit(:,end-1),'*-r','linewidth',2);hold on;
plot(alldkfit(:,6),'*-g','linewidth',2);hold on
% plot(alldkfit(:,6),'*-b','linewidth',2);%%%% should reload data

setgca(18)
xlim([1,size(alldkfit,1)])
subplot(2,1,2);
plot(real(allsigmak2(:,end-1)),'*-r','linewidth',2);hold on
plot(real(allsigmak2(:,6)),'*-g','linewidth',2);hold on
% plot(real(allsigmak2(:,6)),'*-b','linewidth',2);
setgca(18)
xlim([1,size(alldkfit,1)])

lgnd=legend('Exact','update R','constant R');
set(lgnd,'FontSize',16);

%% compare reduced k and full k

load alldkfit;load allsigmak2;load allomegakfit %%% must reload the data
load redind
figure()
subplot(2,1,1);
plot(-alldkfit(redind,end-1),'*-r','linewidth',2);hold on;
plot(alldkfit(redind,6),'*-g','linewidth',2);hold on
% plot(alldkfit(redind,6),'*-b','linewidth',2);%%%% should reload data
xlim([1,length(redind)])
setgca(18)
title('d_k','FontSize',18);
subplot(2,1,2);
plot(real(allsigmak2(redind,end-1)),'*-r','linewidth',2);hold on
plot(real(allsigmak2(redind,6)),'*-g','linewidth',2);hold on
% plot(real(allsigmak2(redind,6)),'*-b','linewidth',2);
setgca(18)
xlim([1,length(redind)])
title('\omega_k','FontSize',18);
lgnd=legend('Exact','full modes','reduced modes');
set(lgnd,'FontSize',16);