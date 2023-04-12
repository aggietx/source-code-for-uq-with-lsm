clear;close all;
cc=[];rms=[];
p=.5;
if p==1
load data_p1_L36_para2d;cc=[cc,ccall];rms=[rms,rmsall];
else
load data_p05_L36_para2d;cc=[cc,ccall];rms=[rms,rmsall];
end

if p==1
load data_p1_L36_para2d_k12;cc=[cc,ccall];rms=[rms,rmsall];
else
load data_p05_L36_para2d_k12;cc=[cc,ccall];rms=[rms,rmsall];
end


figure()
subplot(2,1,1)
plot(1:T,rms(:,1),'g','linewidth',1);hold on
plot(1:T,rms(:,2),'b','linewidth',1);
% lgnd=legend('1D model','2D model','EnKBF 5D','EnKBF perfect');
% set(lgnd,'FontSize',18);
setgca(18)
% xlabel('t','fontsize',24)
xlim([1,T])
title('RMS of velocity','fontsize',18)


subplot(2,1,2)
plot(1:T,cc(:,1),'g','linewidth',1);hold on
plot(1:T,cc(:,2),'b','linewidth',1);
lgnd=legend('Kmax=10','Kmax=12');
set(lgnd,'FontSize',14);
setgca(18)
xlabel('t','fontsize',24)
xlim([1,T])
title('CC of velocity','fontsize',18)

