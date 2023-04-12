clear
% close all
load  data_gb_L12_K5c;
range=0:Nit-1;range1=1:length(exactenekmean);
figure();
subplot(4,5,1)
plot(range,eL2dk(range+1),'-*b','linewidth',1.5);setgca(12)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
% title('Relative L^2 error of d_k, L=12','FontSize',14);
xlim([min(range),max(range)])
ylabel('(I) L=12','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(a) Relative L^2 error of d_k','FontSize',13)

subplot(4,5,6)
plot(range,eL2sigma(range+1),'-*r','linewidth',1.5);setgca(12);xlim([min(range),max(range)])
ylabel('(II) L=12','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(f) Relative L^2 error of \sigma_k','FontSize',13)

subplot(4,5,2)
imagesc([-5,5],[-5,5],dkexact2d);colorbar;temp1=caxis;setgca(12)
title('(b) Exact d_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,3);
imagesc([-5,5],[-5,5],dkfit2d);colorbar;caxis(temp1);setgca(12)
title('(c) Estimated d_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,2+5)
imagesc([-5,5],[-5,5],sigmaexact2d);colorbar;temp1=caxis;setgca(12)
title('(g) Exact \sigma_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,3+5);
imagesc([-5,5],[-5,5],sigmafit2d);colorbar;caxis(temp1);setgca(12)
title('(h) Estimated \sigma_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,4)
imagesc([-5,5],[-5,5],omegafit2d);colorbar;setgca(12);
title('(d) Estimated \omega_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,4+5);
imagesc([-5,5],[-5,5],real(Ffit2d));colorbar;setgca(12)
title('(i) Estimated F_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,[5,5*2])
plot(range1,exactenekmean,'-*k','linewidth',1.5);hold on
plot(range1,estenekmean,'-*g','linewidth',1.5);hold on
plot(range1,sampleenekmean1/79999/2,'-*b','linewidth',1.5);
% subplot(4,5,5*2)
xlim([min(range1),max(range1)]);setgca(12);
title('(e) Energy comparison','FontSize',13)
%%%%%%%%%%%%%
load data_gb_L60_K5c;
subplot(4,5,11)
plot(range,eL2dk(range+1),'-*b','linewidth',1.5);setgca(12);xlim([min(range),max(range)])
ylabel('(III) L=60','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(j) Relative L^2 error of d_k','FontSize',13)

subplot(4,5,16)
plot(range,eL2sigma(range+1),'-*r','linewidth',1.5);setgca(12);xlim([min(range),max(range)])
ylabel('(IV) L=60','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
xlabel('Iteration','FontSize',13)
title('(o) Relative L^2 error of \sigma_k','FontSize',13)

subplot(4,5,2+10)
imagesc([-5,5],[-5,5],dkexact2d);colorbar;temp1=caxis;setgca(12)
title('(k) Exact d_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,3+10);
imagesc([-5,5],[-5,5],dkfit2d);colorbar;caxis(temp1);setgca(12)
title('(l) Estimated d_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,2+5+10)
imagesc([-5,5],[-5,5],sigmaexact2d);colorbar;temp1=caxis;setgca(12)
title('(p) Exact \sigma_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,3+5+10);
imagesc([-5,5],[-5,5],sigmafit2d);colorbar;caxis(temp1);setgca(12)
title('(q) Estimated \sigma_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,4+10)
imagesc([-5,5],[-5,5],omegafit2d);colorbar;setgca(12);
title('(m) Estimated \omega_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

subplot(4,5,4+5+10);
imagesc([-5,5],[-5,5],real(Ffit2d));colorbar;setgca(12)
title('(r) Estimated F_k','FontSize',13)
xlabel('k_x','fontsize',14);
ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)



subplot(4,5,[5+10,5*2+10])
plot(range1,exactenekmean,'-*k','linewidth',1.5);hold on
plot(range1,estenekmean,'-*g','linewidth',1.5);hold on
plot(range1,sampleenekmean1/79999/2,'-*b','linewidth',1.5);
% subplot(4,5,5*2)
xlim([min(range1),max(range1)]);setgca(12);
title('(n) Energy comparison','FontSize',13)
xlabel('|k|','FontSize',13);
lgnd=legend('Exact','Estimated','Reconstructed');
set(lgnd,'FontSize',13);