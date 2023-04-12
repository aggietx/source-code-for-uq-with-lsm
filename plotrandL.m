% close all;
load rcc_gb144b.mat ;  load   rcc_sw144b.mat   
load rcc_gb108b.mat ; load  rcc_gb36b.mat; load    rcc_sw108b.mat  ; load rcc_sw36b.mat    
load rcc_gb144a.mat  ; load rcc_gb72b.mat  ; load  rcc_sw144a.mat  ; load rcc_sw72b.mat 
L=144;
aa=[1:6,12:6:L];ind=[1:6,7:2:29];ind36=[1:6,7:2:11];
ind72=[1:6,7:2:17];ind108=[1:6,7:2:23];
aa(ind);

figure();
subplot(2,1,1);
plot(aa(ind),rcc_gb144a(ind,2),'-*r','linewidth',1.5);hold on
plot(aa(ind36),rcc_gb36b(ind36,2),'-*b','linewidth',1.5);hold on
plot(aa(ind72),rcc_gb72b(ind72,2),'-*g','linewidth',1.5);hold on
plot(aa(ind108),rcc_gb108b(ind108,2),'-*k','linewidth',1.5);hold on
plot(aa(ind),rcc_gb144b(ind,2),'-*c','linewidth',1.5);hold on

setgca(16);
% lgnd=legend('L/L','L/36','L/72','L/108','L/144');
% set(lgnd,'FontSize',14);
title('Mean RMS of velocity (GB)','FontSize',14);
xlim([1,144])
% xlabel('L','FontSize',14);

subplot(2,1,2);
plot([1:6,12:12:144],rcc_sw144a(:,2),'-*r','linewidth',1.5);hold on
plot([1:6,12:12:36],rcc_sw36b(:,2),'-*b','linewidth',1.5);hold on
plot([1:6,12:12:72],rcc_sw72b(:,2),'-*g','linewidth',1.5);hold on
plot([1:6,12:12:108],rcc_sw108b(:,2),'-*k','linewidth',1.5);hold on
plot([1:6,12:12:144],rcc_sw144b(:,2),'-*c','linewidth',1.5);hold on

setgca(16);
lgnd=legend('L/L','L/36','L/72','L/108','L/144');
set(lgnd,'FontSize',14);
title('Mean RMS of velocity (shallow water)','FontSize',14);
xlim([1,144])
xlabel('L','FontSize',14);
return
%%%%%%%%% local data%%%%%%%%%%%%
err=[ 14.00 2.40  0.30 %%% 12
 6.52  1.11  0.54; %%%18
 6.15  1.03  0.61  %% 24
 3.49  0.58  0.82%%%% 30
 2.78  0.47  0.88; %%%34
 2.42,0.41,0.91];%%36
ll=[12,18,24,30,34,36];
figure();
% subplot(2,1,1);
plot(ll,err(:,2),'-*r','linewidth',1.5);hold on

setgca(16);
% lgnd=legend('L/L','L/36','L/72','L/108','L/144');
% set(lgnd,'FontSize',14);
title('Mean RMS of velocity (GB localization)','FontSize',14);
xlim([12,36])
xlabel('L','FontSize',14);

err=[ 328.27 91.44  0.02 %%%% r=2
    6.66 1.12  0.54  %%% r=2.5
    5.96 0.99  0.62 %%% r=3
   2.94  0.49  0.87                   %%% r=3.5
   2.55, 0.43  0.90];%r=4
ll=[2,2.5,3,3.5,4];
figure();
% subplot(2,1,1);
plot(ll(2:end),err(2:end,2),'-*r','linewidth',1.5);hold on

setgca(16);
% lgnd=legend('L/L','L/36','L/72','L/108','L/144');
% set(lgnd,'FontSize',14);
title('Mean RMS of velocity (GB localization)','FontSize',14);
xlim([2,4])
xlabel('r','FontSize',14);
