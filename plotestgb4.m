clear ;close all;
%%% 2d parameter comparison
% eL2dk; eL2omega;eL2F;eL2sigma;exactpara2d;estpara2d;frdata;savedtraj
% savedtraj=[u_hat(:,gap1:gap1:end);
%     exact_gamma_mean_trace_smoother(:,1:gap1:end);
%     gamma_mean_trace_smoother(:,1:gap1:end)];
% % set(gca,'XTick',[12,24,48,72,96],'FontSize',20)

% load data_gb_L12_K5.mat;
load data_gb_L36_K5.mat;
% a=exactpara2d(:,:,4);
% b=estpara2d(:,:,4);
% plot2line(a(:),b(:))
ix=k_0;iy=k_0;ind=getind(ix,iy,kk);
savedtraj=real(savedtraj);
subplot(2,1,1)
plot(1:T,savedtraj(ind,:),'r','linewidth',1.5);hold on;
plot(1:T,savedtraj(ind+dimuhat,:),'b','linewidth',1.5);hold on;
plot(1:T,savedtraj(ind+dimuhat*2,:),'g','linewidth',1.5);
setgca(16);

title(['GB mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',14);

ix=Kmax;iy=Kmax;ind=getind(ix,iy,kk);
subplot(2,1,2)
plot(1:T,savedtraj(ind,:),'r','linewidth',1.5);hold on;
plot(1:T,savedtraj(ind+dimuhat,:),'b','linewidth',1.5);hold on;
plot(1:T,savedtraj(ind+dimuhat*2,:),'g','linewidth',1.5);
setgca(16);
lgnd=legend('Truth','Exact paramter','Estimated paramter');
set(lgnd,'FontSize',12);
xlabel('t','FontSize',16);
title(['GB mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',14);
