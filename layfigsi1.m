%%% figure in supplementary 
%%%% plot recovered energy
% load data2d_layer1_L6_p2_K7
clear;
load enereg2;
load enereg1;
close all

figure();
 load  data2d_layer1_L6_p2.mat 
 
 F1=real(allFfit(:,Nit));F2=imag(allFfit(:,Nit));%%%% use estimated parameter to compute the energy
 estene=1/2*((F2.^2+ F1.^2)./alldkfit(:,Nit) +allsigmak2(:,Nit).^2*2./(2*alldkfit(:,Nit)) );
%   estene1=1/2*(allsigmak2(:,Nit).^2*2./(2*alldkfit(:,Nit)) );
% estene=getene(estuhat);
 ene2d=plotest2d(Kmax,kk,estene);
 
 eneexact2d=10^(-3.8)*ones(13,13);
eneexact2d(7,7:end)=enereg1;
eneexact2d(7,1:7)=enereg1(end:-1:1);
subplot(2,2,1);
  imagesc([-5,5],[-5,5],log10(eneexact2d));colorbar;temp1=caxis;
 setgca(16)
title('(a) Regime I: Nearly Gaussian','FontSize',12)

 
  subplot(2,2,3)
  imagesc([-5,5],[-5,5],log10(ene2d));colorbar;caxis(temp1);
%   temp1=caxis;
  setgca(16)
% xlabel('k_x','fontsize',14);
% ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
% plot(0:6,enereg1,'*-k','linewidth',2);hold on
% plot(0:6,ene2d(7,7:end),'*-b','linewidth',2);hold on
% setgca(16);
% xlabel('k_x','fontsize',16');
% xlim([0,6])

   load  data2d_layer2_L6_p2.mat 
 
 F1=real(allFfit(:,Nit));F2=imag(allFfit(:,Nit));%%%% use estimated parameter to compute the energy
 estene=1/2*((F2.^2+ F1.^2)./alldkfit(:,Nit) +allsigmak2(:,Nit).^2*2./(2*alldkfit(:,Nit)) );
%  estene=getene(estuhat);
 ene2d2=plotest2d(Kmax,kk,estene);

 eneexact2d=10^(-2.5)*ones(13,13);
eneexact2d(7,7:end)=enereg2;
eneexact2d(7,1:7)=enereg2(end:-1:1);
subplot(2,2,2);
  imagesc([-5,5],[-5,5],log10(eneexact2d));colorbar;caxis([-2.6,0.92]);temp=caxis;
 setgca(16)
 title('(b) Regime II: Strongly non-Gaussian','FontSize',12)

    subplot(2,2,4)
  imagesc([-5,5],[-5,5],log10(ene2d2));colorbar;caxis(temp);
%   temp1=caxis;
  setgca(16)
% xlabel('k_x','fontsize',14);
% ylabel('k_y','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
% enereg1(1)=enereg1(1)*2;
% enereg1(2:end)=enereg1(2:end)*4;



  
% subplot(2,2,4);
% plot(0:6,enereg2,'*-k','linewidth',2);hold on
% plot(0:6,ene2d2(7,7:end),'*-b','linewidth',2);hold on
% setgca(16);
% xlabel('k_x','fontsize',16');
% xlim([0,6])
% 
% lgnd=legend('Truth','Recovered spectrums w/ est param');
% set(lgnd,'FontSize',12);