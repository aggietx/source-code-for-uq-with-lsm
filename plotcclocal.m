load allccfilgb;load allrmsfilgb
load allcc_gbr3;load allrms_gbr3
load allcc_gbr25;load allrms_gbr25
dt=1/200;T=50;
figure()
subplot(2,1,1)
plot(dt:dt:T,allrmsfilgb,'b','linewidth',1);hold on
plot(dt:dt:T,allrms_gbr3,'g','linewidth',1);hold on
plot(dt:dt:T,allrms_gbr25,'r','linewidth',1);
% lgnd=legend('1D model','2D model','EnKBF 5D','EnKBF perfect');
% set(lgnd,'FontSize',18);
setgca(18)
% xlabel('t','fontsize',24)
xlim([1,T])
title('RMS of velocity','fontsize',18)


subplot(2,1,2)
plot(dt:dt:T,allccfilgb,'b','linewidth',1);hold on
plot(dt:dt:T,allcc_gbr3,'g','linewidth',1);hold on
plot(dt:dt:T,allcc_gbr25,'r','linewidth',1);
lgnd=legend('No localization','r=3 (L=25.8)','r=2.5 (L=17.9)');
set(lgnd,'FontSize',14);
setgca(18)
xlabel('t','fontsize',24)
xlim([1,T])
title('CC of velocity','fontsize',18)
fprintf('mean of cc %2.2f %2.2f %2.2f\n',mean(allccfilgb),mean(allcc_gbr3),mean(allcc_gbr25));