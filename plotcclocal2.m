load allccfilsw;load allrmsfilsw
load allcc_swr32;load allrms_swr32
load allcc_swr3;load allrms_swr3
dt=1/200;T=20;
figure()
subplot(2,1,1)
plot(dt:dt:T,allrmsfilsw,'b','linewidth',1);hold on
plot(dt:dt:T,allrms_swr32,'g','linewidth',1);hold on
plot(dt:dt:T,allrms_swr3,'r','linewidth',1);
% lgnd=legend('1D model','2D model','EnKBF 5D','EnKBF perfect');
% set(lgnd,'FontSize',18);
setgca(18)
% xlabel('t','fontsize',24)
xlim([1,T])
title('RMS of velocity','fontsize',18)


subplot(2,1,2)
plot(dt:dt:T,allccfilsw,'b','linewidth',1);hold on
plot(dt:dt:T,allcc_swr32,'g','linewidth',1);hold on
plot(dt:dt:T,allcc_swr3,'r','linewidth',1);
lgnd=legend('No localization','r=3.2 (L=29.2)','r=3 (L=17.9)');
set(lgnd,'FontSize',14);
setgca(18)
xlabel('t','fontsize',24)
xlim([1,T])
title('CC of velocity','fontsize',18)
fprintf('mean of cc %2.2f %2.2f %2.2f\n',mean(allccfilsw),mean(allcc_swr32),mean(allcc_swr3));