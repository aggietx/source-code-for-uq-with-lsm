clear ;close all;
load data_sw_L36_K3.mat;id=1;rmst(id,:)=errest(1,1:Nit);rms(id,:)=errest(2,1:Nit);cc(id,:)=errest(3,1:Nit);
errsm=[mean(rmstimeexactsm(3:end));mean(rmsexactsm);mean(ccexactsm)];errsmL(:,id)=errsm;
id=1;freermst(id)=errest(1,end);freerms(id)=errest(2,end);freecc(id)=errest(3,end);

load data_sw_L72_K3.mat;id=2;rmst(id,:)=errest(1,1:Nit);rms(id,:)=errest(2,1:Nit);cc(id,:)=errest(3,1:Nit);
errsm=[mean(rmstimeexactsm(3:end));mean(rmsexactsm);mean(ccexactsm)];errsmL(:,id)=errsm;
freermst(id)=errest(1,end);freerms(id)=errest(2,end);freecc(id)=errest(3,end);

load data_sw_L108_K3.mat;id=3;rmst(id,:)=errest(1,1:Nit);rms(id,:)=errest(2,1:Nit);cc(id,:)=errest(3,1:Nit);
errsm=[mean(rmstimeexactsm(3:end));mean(rmsexactsm);mean(ccexactsm)];errsmL(:,id)=errsm;
freermst(id)=errest(1,end);freerms(id)=errest(2,end);freecc(id)=errest(3,end);

load data_sw_L144_K3.mat;id=4;rmst(id,:)=errest(1,1:Nit);rms(id,:)=errest(2,1:Nit);cc(id,:)=errest(3,1:Nit);
errsm=[mean(rmstimeexactsm(3:end));mean(rmsexactsm);mean(ccexactsm)];errsmL(:,id)=errsm;
freermst(id)=errest(1,end);freerms(id)=errest(2,end);freecc(id)=errest(3,end);

load data_sw_L180_K3.mat;id=5;rmst(id,:)=errest(1,1:Nit);rms(id,:)=errest(2,1:Nit);cc(id,:)=errest(3,1:Nit);
errsm=[mean(rmstimeexactsm(3:end));mean(rmsexactsm);mean(ccexactsm)];errsmL(:,id)=errsm;
freermst(id)=errest(1,end);freerms(id)=errest(2,end);freecc(id)=errest(3,end);

length=size(cc,2);range=0:length-1;range2=1:length;
figure();
subplot(3,1,1)
plot(range,rmst(1,range2),'-*r','linewidth',1.5);hold on;
plot(range,rmst(2,range2),'-*b','linewidth',1.5);hold on;
plot(range,rmst(3,range2),'-*green','linewidth',1.5);hold on;
plot(range,rmst(4,range2),'-*k','linewidth',1.5);hold on;
plot(range,rmst(5,range2),'-*c','linewidth',1.5);hold on;
% plot(range,rms(6,range2),'-*y','linewidth',1.5);hold on;
setgca(16)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
title('Mean RMS of velocity over time','FontSize',14);
xlim([min(range),max(range)])



subplot(3,1,2)
plot(range,rms(1,range2),'-*r','linewidth',1.5);hold on;
plot(range,rms(2,range2),'-*b','linewidth',1.5);hold on;
plot(range,rms(3,range2),'-*green','linewidth',1.5);hold on;
plot(range,rms(4,range2),'-*k','linewidth',1.5);hold on;
plot(range,rms(5,range2),'-*c','linewidth',1.5);hold on;
% plot(range,rms(6,range2),'-*y','linewidth',1.5);hold on;
setgca(16)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
title('Mean RMS of velocity over fixed point','FontSize',14);
xlim([min(range),max(range)])

subplot(3,1,3)
plot(range,cc(1,range2),'-*r','linewidth',1.5);hold on;
plot(range,cc(2,range2),'-*b','linewidth',1.5);hold on;
plot(range,cc(3,range2),'-*green','linewidth',1.5);hold on;
plot(range,cc(4,range2),'-*k','linewidth',1.5);hold on;
plot(range,cc(5,range2),'-*c','linewidth',1.5);hold on;
% plot(range,cc(6,range2),'-*y','linewidth',1.5);hold on;
setgca(16)
lgnd=legend('L=36','L=72','L=108','L=144','L=180');
set(lgnd,'FontSize',12);
title('Mean CC of velocity over fixed point','FontSize',14);
xlim([min(range),max(range)])
xlabel('Iteration','FontSize',16);


range3=36:36:180;
% range3=1:5;
figure();
subplot(3,1,1)
plot(range3,errsmL(1,:),'-*r','linewidth',1.5);hold on;
plot(range3,rmst(:,end),'-*b','linewidth',1.5);hold on;
plot(range3,freermst,'-*g','linewidth',1.5);
setgca(16)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
title('Mean RMS of velocity over time','FontSize',14);
xlim([min(range3),max(range3)])



subplot(3,1,2);
id=2;
plot(range3,errsmL(id,:),'-*r','linewidth',1.5);hold on;
plot(range3,rms(:,end),'-*b','linewidth',1.5);hold on;
plot(range3,freerms,'-*g','linewidth',1.5);
setgca(16)
% lgnd=legend('Exact','Estimation');
% set(lgnd,'FontSize',14);
title('Mean RMS of velocity over fixed point','FontSize',14);
xlim([min(range3),max(range3)])

subplot(3,1,3)
id=3;
plot(range3,errsmL(id,:),'-*r','linewidth',1.5);hold on;
plot(range3,cc(:,end),'-*b','linewidth',1.5);hold on;
plot(range3,freecc,'-*g','linewidth',1.5);
setgca(16)
lgnd=legend('Exact parameter','estimated paramter','free run with fixed noise');
set(lgnd,'FontSize',12);
title('Mean CC of velocity over fixed point','FontSize',14);
xlim([min(range3),max(range3)])
xlabel('Iteration','FontSize',16);