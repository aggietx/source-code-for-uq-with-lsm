%%%% filter gb data 
clear;close all;
load datafil_L72_sw_K3_case1;velerrcase1(:,1)=errfil;
load datafil_L72_sw_K3_case2;velerrcase2(:,1)=errfil;
load datafil_L72_sw_K3_case3;velerrcase3(:,1)=errfil;

load datafil_L108_sw_K3_case1;velerrcase1(:,2)=errfil;
load datafil_L108_sw_K3_case2;velerrcase2(:,2)=errfil;
load datafil_L108_sw_K3_case3;velerrcase3(:,2)=errfil;

load datafil_L144_sw_K3_case1;velerrcase1(:,3)=errfil;
load datafil_L144_sw_K3_case2;velerrcase2(:,3)=errfil;
load datafil_L144_sw_K3_case3;velerrcase3(:,3)=errfil;

load datafil_L180_sw_K3_case1;velerrcase1(:,4)=errfil;
load datafil_L180_sw_K3_case2;velerrcase2(:,4)=errfil;
load datafil_L180_sw_K3_case3;velerrcase3(:,4)=errfil;


nL=12*6:12*3:12*15;
figure()
subplot(3,1,1)
plot(nL,velerrcase1(1,:),'-*g','linewidth',1.5);hold on
plot(nL,velerrcase2(1,:),'-*b','linewidth',1.5);hold on
plot(nL,velerrcase3(1,:),'-*r','linewidth',1.5);hold on
% lgnd=legend('1D model','2D model','EnKBF Reduced','EnKBF perfect','Exact Reduced');
% set(lgnd,'fontsize',16;
setgca(14)
% xlabel('L','fontsize',14)
xlim([min(nL),max(nL)])
title('Mean RMS of velocity over time','fontsize',14)
% lgnd=legend('Regular R','Diagonal R','Constant R');
% set(lgnd,'FontSize',14);

subplot(3,1,2)
plot(nL,velerrcase1(2,:),'-*g','linewidth',1.5);hold on
plot(nL,velerrcase2(2,:),'-*b','linewidth',1.5);hold on
plot(nL,velerrcase3(2,:),'-*r','linewidth',1.5);hold on
% lgnd=legend('1D model','2D model','EnKBF Reduced','EnKBF perfect','Exact Reduced');
% set(lgnd,'fontsize',16;
setgca(14)
% xlabel('L','fontsize',14)
xlim([min(nL),max(nL)])
title('Mean RMS of velocity over fixed point','fontsize',14)
% lgnd=legend('Regular R','Diagonal R','Constant R');
% set(lgnd,'FontSize',14);

subplot(3,1,3)
plot(nL,velerrcase1(3,:),'-*g','linewidth',1.5);hold on
plot(nL,velerrcase2(3,:),'-*b','linewidth',1.5);hold on
plot(nL,velerrcase3(3,:),'-*r','linewidth',1.5);hold on
% lgnd=legend('1D model','2D model','EnKBF Reduced','EnKBF perfect','Exact Reduced');
% set(lgnd,'fontsize',16;
setgca(14)
xlabel('L','fontsize',14)
xlim([min(nL),max(nL)])
title('Mean CCof velocity over fixed point','fontsize',14)
lgnd=legend('Regular R','Diagonal R','Constant R');
set(lgnd,'FontSize',14);