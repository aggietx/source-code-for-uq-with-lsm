%%%% filter gb data 
clear;close all;
load datafil_L12_gb_K5_case1;velerrcase1(:,1)=errfil;
load datafil_L12_gb_K5_case2;velerrcase2(:,1)=errfil;
load datafil_L12_gb_K5_case3;velerrcase3(:,1)=errfil;

load datafil_L24_gb_K5_case1;velerrcase1(:,2)=errfil;
load datafil_L24_gb_K5_case2;velerrcase2(:,2)=errfil;
load datafil_L24_gb_K5_case3;velerrcase3(:,2)=errfil;

load datafil_L36_gb_K5_case1;velerrcase1(:,3)=errfil;
load datafil_L36_gb_K5_case2;velerrcase2(:,3)=errfil;
load datafil_L36_gb_K5_case3;velerrcase3(:,3)=errfil;

load datafil_L48_gb_K5_case1;velerrcase1(:,4)=errfil;
load datafil_L48_gb_K5_case2;velerrcase2(:,4)=errfil;
load datafil_L48_gb_K5_case3;velerrcase3(:,4)=errfil;

load datafil_L60_gb_K5_case1;velerrcase1(:,5)=errfil;
load datafil_L60_gb_K5_case2;velerrcase2(:,5)=errfil;
load datafil_L60_gb_K5_case3;velerrcase3(:,5)=errfil;

load datafil_L72_gb_K5_case1;velerrcase1(:,6)=errfil;
load datafil_L72_gb_K5_case2;velerrcase2(:,6)=errfil;
load datafil_L72_gb_K5_case3;velerrcase3(:,6)=errfil;
nL=12:12:72;
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