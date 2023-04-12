%%%% dispersion and signal, and energy
%%%% sigest sigsm, disest, dissm
%%%% relative entropy
clear ;close all;st=5;
load data_gb_L12_K5.mat;id=1;sig(id,:)=[mean(sigest(st:end)),mean(sigsm(st:end))];dis(id,:)=[mean(disest(st:end)),mean(dissm(st:end))];
load data_gb_L24_K5.mat;id=2;sig(id,:)=[mean(sigest(st:end)),mean(sigsm(st:end))];dis(id,:)=[mean(disest(st:end)),mean(dissm(st:end))];
load data_gb_L36_K5.mat;id=3;sig(id,:)=[mean(sigest(st:end)),mean(sigsm(st:end))];dis(id,:)=[mean(disest(st:end)),mean(dissm(st:end))];
load data_gb_L48_K5.mat;id=4;sig(id,:)=[mean(sigest(st:end)),mean(sigsm(st:end))];dis(id,:)=[mean(disest(st:end)),mean(dissm(st:end))];
load data_gb_L60_K5.mat;id=5;sig(id,:)=[mean(sigest(st:end)),mean(sigsm(st:end))];dis(id,:)=[mean(disest(st:end)),mean(dissm(st:end))];
range3=12:12:60;
% range3=1:5;
figure();
subplot(2,1,1)
plot(range3,dis(:,2),'-*r','linewidth',1.5);hold on;
plot(range3,dis(:,1),'-*b','linewidth',1.5);hold on;
setgca(16)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
title('Dispersion','FontSize',14);
xlim([min(range3),max(range3)])
setgca(16)
subplot(2,1,2)
plot(range3,sig(:,2),'-*r','linewidth',1.5);hold on;
plot(range3,sig(:,1),'-*b','linewidth',1.5);hold on;
lgnd=legend('Exact parameter','estimated paramter');

title('Signal','FontSize',18);
xlim([min(range3),max(range3)])
setgca(16)

xlabel('L','FontSize',16);