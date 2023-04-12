clear ;close all;
%%% 2d parameter comparison
% eL2dk; eL2omega;eL2F;eL2sigma;exactpara2d;estpara2d
load data_gb_L36_K5.mat;
% a=exactpara2d(:,:,4);
% b=estpara2d(:,:,4);
% plot2line(a(:),b(:))

figure();
subplot(2,2,1);
imagesc(exactpara2d(:,:,1));colorbar;setgca(14)
title('Exact d_k','FontSize',14);

subplot(2,2,2);
imagesc(estpara2d(:,:,1));colorbar;setgca(14)
title('Estimated d_k','FontSize',14);

subplot(2,2,3);
imagesc(exactpara2d(:,:,4));colorbar;setgca(14)
title('Estimated d_k','FontSize',14);
subplot(2,2,4);
imagesc(estpara2d(:,:,4));colorbar;setgca(14)
title('Estimated \omega_k','FontSize',14);
