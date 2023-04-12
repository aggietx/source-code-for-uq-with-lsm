clear;
% close all;
load xexact1;load yexact1.mat
load dataswcase1_L30;
S=1.5;figure();
rmstimefil;
rmstimefilgb;
rmstimefilgr;
good=1;
if good
time1=find(rmstimefilgr==min(rmstimefilgr));
min(rmstimefilgr)
else
  time1=find(rmstimefilgr==max(rmstimefilgr)); 
  max(rmstimefilgr)
end
data=(savedtrajsw);
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,25,data(1:dimuhat,time1));
subplot(3,3,1)
% imagesc(sqrt(ux.^2+uy.^2));lightcolor; temp1=caxis;hold on
contourf(xx,yy,sqrt(ux.^2+uy.^2),'LineStyle','none');lightcolor;temp1=caxis;hold on
quiver(xx,yy,ux,uy,S, 'linewidth',1,'color','k');
hold on;scatter(xexact1(1:30,time1),yexact1(1:30,time1),'ro','filled');
setgca(13);xlim([0,L1]);ylim([0,L1])
ylabel('(I) Truth,','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(a) SW','FontSize',13);

subplot(3,3,2)
data2=data(1:dimuhat,time1);data2(1+dimuhat0:end,:)=0;
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,25,data2);
% imagesc(sqrt(ux.^2+uy.^2));lightcolor; temp2=caxis;hold on
contourf(xx,yy,sqrt(ux.^2+uy.^2),'LineStyle','none');lightcolor;temp2=caxis;hold on
quiver(xx,yy,ux,uy,S, 'linewidth',1,'color','k');
hold on;scatter(xexact1(1:30,time1),yexact1(1:30,time1),'ro','filled');
setgca(13);xlim([0,L1]);ylim([0,L1])
% ylabel('(I) Truth,','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(b) GB of SW','FontSize',13);

subplot(3,3,3)
data2=data(1:dimuhat,time1);data2(1:dimuhat0,:)=0;
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,25,data2);
% imagesc(sqrt(ux.^2+uy.^2));lightcolor; temp3=caxis;hold on
contourf(xx,yy,sqrt(ux.^2+uy.^2),'LineStyle','none');lightcolor; temp3=caxis;hold on
quiver(xx,yy,ux,uy,S, 'linewidth',1,'color','k');
hold on;scatter(xexact1(1:30,time1),yexact1(1:30,time1),'ro','filled');
setgca(13);xlim([0,L1]);ylim([0,L1])
% ylabel('(I) Truth,','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(c) GR of SW','FontSize',13);

%%%%%%
subplot(3,3,4);
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,25,data(1+dimuhat:end,time1));
% imagesc(sqrt(ux.^2+uy.^2));lightcolor; caxis(temp1);hold on
contourf(xx,yy,sqrt(ux.^2+uy.^2),'LineStyle','none'); lightcolor;caxis(temp1);hold on
quiver(xx,yy,ux,uy,S, 'linewidth',1,'color','k');
% hold on;scatter(xexact1(1:30,time1),yexact1(1:30,time1),'ro','filled');
setgca(13);xlim([0,L1]);ylim([0,L1])
ylabel('(II) Full R,','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)

subplot(3,3,5)
data2=data(1+dimuhat:end,time1);data2(1+dimuhat0:end,:)=0;
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,25,data2);
% imagesc(sqrt(ux.^2+uy.^2));lightcolor; caxis(temp2);hold on
contourf(xx,yy,sqrt(ux.^2+uy.^2),'LineStyle','none');lightcolor;caxis(temp2);hold on
quiver(xx,yy,ux,uy,S, 'linewidth',1,'color','k');
% hold on;scatter(xexact1(1:30,time1),yexact1(1:30,time1),'ro','filled');
setgca(13);xlim([0,L1]);ylim([0,L1])
% ylabel('(I) Truth,','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('(b) GB of SW','FontSize',13);

subplot(3,3,6)
data2=data(1+dimuhat:end,time1);data2(1:dimuhat0,:)=0;
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,25,data2);
% imagesc(sqrt(ux.^2+uy.^2));lightcolor; caxis(temp3);hold on
contourf(xx,yy,sqrt(ux.^2+uy.^2),'LineStyle','none');lightcolor; caxis(temp3);hold on
quiver(xx,yy,ux,uy,S, 'linewidth',1,'color','k');
% hold on;scatter(xexact1(1:30,time1),yexact1(1:30,time1),'ro','filled');
setgca(13);xlim([0,L1]);ylim([0,L1])
% ylabel('(I) Truth,','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('(b) GR of SW','FontSize',13);


%%%%%%%%%
load dataswcase3_L30;data=(savedtrajsw);
subplot(3,3,7)
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,25,data(1+dimuhat:end,time1));
% imagesc(sqrt(ux.^2+uy.^2));lightcolor; caxis(temp1);colorbar('horzontial');    hold on
contourf(xx,yy,sqrt(ux.^2+uy.^2),'LineStyle','none');lightcolor; caxis(temp1);colorbar('horzontial');    hold on
quiver(xx,yy,ux,uy,S, 'linewidth',1,'color','k');
% hold on;scatter(xexact1(1:30,time1),yexact1(1:30,time1),'ro','filled');
setgca(13);xlim([0,L1]);ylim([0,L1])
ylabel('(III) Constant R,','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)

subplot(3,3,8)
data2=data(1+dimuhat:end,time1);data2(1+dimuhat0:end,:)=0;
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,25,data2);
% imagesc(sqrt(ux.^2+uy.^2));lightcolor; caxis(temp2);colorbar('horzontial');    hold on
contourf(xx,yy,sqrt(ux.^2+uy.^2),'LineStyle','none'); lightcolor;caxis(temp2);colorbar('horzontial');    hold on
quiver(xx,yy,ux,uy,S, 'linewidth',1,'color','k');
% hold on;scatter(xexact1(1:30,time1),yexact1(1:30,time1),'ro','filled');
setgca(13);xlim([0,L1]);ylim([0,L1])
% ylabel('(I) Truth,','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('(b) GB of SW','FontSize',13);

subplot(3,3,9)
data2=data(1+dimuhat:end,time1);data2(1:dimuhat0,:)=0;
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,25,data2);
% imagesc(sqrt(ux.^2+uy.^2));lightcolor; caxis(temp3);colorbar('horzontial');    hold on
contourf(xx,yy,sqrt(ux.^2+uy.^2),'LineStyle','none');lightcolor;caxis(temp3);colorbar('horzontial');    hold on
quiver(xx,yy,ux,uy,S, 'linewidth',1,'color','k');
% hold on;scatter(xexact1(1:30,time1),yexact1(1:30,time1),'ro','filled');
setgca(13);xlim([0,L1]);ylim([0,L1])
% ylabel('(I) Truth,','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('(b) GR of SW','FontSize',13);


% load dataswcase3_L30

