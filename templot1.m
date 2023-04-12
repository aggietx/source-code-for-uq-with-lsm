%%%% plot 1d data
clear;close all;
load ccall1d
load ccall2d
load estuhat ;load uhate;load estuhat2;
figure()
plot(ccall2d,'-*r','linewidth',1);
hold on;
plot(ccall1d,'-*b','linewidth',1);
setgca(18)
lgnd=legend('2d','1d');
set(lgnd,'FontSize',16);
title('cc of velocity','FontSize',16);
% corrcoef(ccall2d,ccall1d)
xlabel('t','FontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max(abs(imag(estuhat(1,:))))
estuhat=real(estuhat);
% plot2line(uhate(1,:),real(estuhat(1,:)));
figure()
plot(uhate(1,:),'-*r','linewidth',1);
hold on;
plot(estuhat(1,:),'-*b','linewidth',1);
setgca(18)
lgnd=legend('Exact','Estimation');
set(lgnd,'FontSize',16);
title('u','FontSize',16);
rmscc(real(estuhat(1,:)),real(uhate(1,:)),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(real(uhate(2,:)),'-*r','linewidth',1);
hold on;
plot(real(estuhat(2,:)),'-*b','linewidth',1);
setgca(18)
lgnd=legend('Exact','Estimation');
set(lgnd,'FontSize',16);
title('\psi_1','FontSize',16);

rmscc(real(estuhat(2,:)),real(uhate(2,:)),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(real(uhate(4,:)),'-*r','linewidth',1);
hold on;
plot(real(estuhat(4,:)),'-*b','linewidth',1);
setgca(18)
lgnd=legend('Exact','Estimation');
set(lgnd,'FontSize',16);
title('\psi_2','FontSize',16);

rmscc(real(estuhat(4,:)),real(uhate(4,:)),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
i=3;
plot(real(uhate(2*i,:)),'-*r','linewidth',1);
hold on;
plot(real(estuhat(2*i,:)),'-*b','linewidth',1);
setgca(18)
lgnd=legend('Exact','Estimation');
set(lgnd,'FontSize',16);
% title('\psi_2','FontSize',16);

rmscc(real(estuhat(2*i,:)),real(uhate(2*i,:)),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%reduced model
figure()
plot(uhate(1,:),'-*r','linewidth',1);
hold on;
plot(estuhat2(1,:),'-*b','linewidth',1);
setgca(18)
lgnd=legend('Exact','Estimation');
set(lgnd,'FontSize',16);
title('u','FontSize',16);
rmscc(real(estuhat2(1,:)),real(uhate(1,:)),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
i=1;
plot(real(uhate(2*i,:)),'-*r','linewidth',1);
hold on;
plot(real(estuhat2(2*i,:)),'-*b','linewidth',1);
setgca(18)
lgnd=legend('Exact','Estimation');
set(lgnd,'FontSize',16);
title('\psi_1','FontSize',16);

