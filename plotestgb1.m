clear ;close all;
% eL2dk; eL2omega;eL2F;eL2sigma;exactpara2d;estpara2d frdata
load data_gb_L12_K5.mat;
% load data_gb_L36_K5.mat;
range=0:Nit-1;
figure();
subplot(2,2,1)
plot(range,eL2dk(range+1),'-*r','linewidth',1.5);hold on;
setgca(16)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
title('Relative L^2 error of d_k, L=12','FontSize',14);
xlim([min(range),max(range)])


subplot(2,2,2)
plot(range,eL2sigma(range+1),'-*r','linewidth',1.5);hold on;
setgca(16)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
title('Relative L^2 error of \omega_k, L=12','FontSize',14);
xlim([min(range),max(range)])

% xlabel('Iteration','FontSize',16);


 load data_gb_L36_K5.mat;
subplot(2,2,3)
plot(range,eL2dk(range+1),'-*r','linewidth',1.5);hold on;
setgca(16)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
title('Relative L^2 error of d_k, L=36','FontSize',14);
xlim([min(range),max(range)])
xlabel('Iteration','FontSize',16);

subplot(2,2,4)
plot(range,eL2sigma(range+1),'-*r','linewidth',1.5);hold on;
setgca(16)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
title('Relative L^2 error of \omega_k, L=36','FontSize',14);
xlim([min(range),max(range)])

xlabel('Iteration','FontSize',16);