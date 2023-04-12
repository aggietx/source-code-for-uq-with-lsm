

clear
% close all
load  data_gb_L12_K5;
figure();range=1:T;
% timeplot1=100;timeplot2=200;
data=savedtraj(getind(0,1,kk),:);data=real(data);timeplot1=find(data==max(data));
data=savedtraj(getind(2,2,kk),:);data=real(data);timeplot2=find(data==min(data));

subplot(4,3,1);
ix=0;iy=1;ind=getind(ix,iy,kk);
plot(range,real(savedtraj(ind,:)),'k','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat,:)),'r','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat*2,:)),'b','linewidth',1);hold on;
plot(timeplot1,real(savedtraj(ind,timeplot1)),'*r','linewidth',4);
% xlim([min(range),max(range)])
xlim([1,400])
setgca(14)
title(['(a) Mode (', num2str(ix),' , ', num2str(iy), ')'],'FontSize',13);
ylabel('(I) L=12','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)


subplot(4,3,2);
ix=0;iy=2;ind=getind(ix,iy,kk);
plot(range,real(savedtraj(ind,:)),'k','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat,:)),'r','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat*2,:)),'b','linewidth',1);hold on;
% xlim([min(range),max(range)])
xlim([1,400])
setgca(14);title(['(b) Mode (', num2str(ix),' , ', num2str(iy), ')'],'FontSize',13);



subplot(4,3,3);
ix=2;iy=2;ind=getind(ix,iy,kk);
plot(range,real(savedtraj(ind,:)),'k','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat,:)),'r','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat*2,:)),'b','linewidth',1);hold on;
plot(timeplot2,real(savedtraj(ind,timeplot2)),'*r','linewidth',4);
% xlim([min(range),max(range)])
xlim([1,400])
setgca(14);title(['(c) Mode (', num2str(ix),' , ', num2str(iy), ')'],'FontSize',13);


subplot(4,3,1+3*2)
data1=savedtraj(1:dimuhat,timeplot1);
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,Dim_Grid,data1);
l = streamslice(xx,yy,ux,uy);set(l,'LineWidth',1);xlabel('x','fontsize',16);%ylabel('yy','fontsize',16);set(get(gca,'YLabel'),'Rotation',0)
axis tight
set(l,'LineWidth',1)
set(l,'Color','k');
setgca(14);
ylabel(['(III) t=',num2str(timeplot1)],'FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(g) Truth','fontsize',16)

subplot(4,3,2+3*2)
data1=savedtraj(1+2*dimuhat:3*dimuhat,timeplot1);
[ux,uy]=computeGBvel(rk,kk,L1,Dim_Grid,data1);
l = streamslice(xx,yy,ux,uy);
axis tight
set(l,'LineWidth',1);%xlabel('x','fontsize',16);ylabel('y','fontsize',16);set(get(gca,'YLabel'),'Rotation',0)
set(l,'Color','b');
setgca(14);
title('(h) Estimated with L=12','fontsize',16)

subplot(4,3,2+3*3)
data1=savedtraj(1+2*dimuhat:3*dimuhat,timeplot2);
[ux,uy]=computeGBvel(rk,kk,L1,Dim_Grid,data1);
l = streamslice(xx,yy,ux,uy);%xlabel('x','fontsize',16);ylabel('y','fontsize',16);set(get(gca,'YLabel'),'Rotation',0)
axis tight
set(l,'LineWidth',1);%xlabel('x','fontsize',16);ylabel('y','fontsize',16);set(get(gca,'YLabel'),'Rotation',0)
set(l,'Color','b');
setgca(14);
title('(k) Estimated with L=12','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load  data_gb_L60_K5;ix=0;iy=1;ind=getind(ix,iy,kk);
subplot(4,3,3*1+1);
plot(range,real(savedtraj(ind,:)),'k','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat,:)),'r','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat*2,:)),'b','linewidth',1);
xlabel('t','fontsize',16)

% xlim([min(range),max(range)])
xlim([1,400])
setgca(14);title(['(d) Mode (', num2str(ix),' , ', num2str(iy), ')'],'FontSize',13);
ylabel('(II) L=60','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
% title('Truth','fontsize',16)

subplot(4,3,2+3*1);
ix=0;iy=2;ind=getind(ix,iy,kk);
plot(range,real(savedtraj(ind,:)),'k','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat,:)),'r','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat*2,:)),'b','linewidth',1);hold on;
% xlim([min(range),max(range)])
xlim([1,400])
setgca(14);title(['(e) Mode (', num2str(ix),' , ', num2str(iy), ')'],'FontSize',13);
xlabel('t','fontsize',16)




subplot(4,3,3+3*1);
ix=2;iy=2;ind=getind(ix,iy,kk);
plot(range,real(savedtraj(ind,:)),'k','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat,:)),'r','linewidth',1);hold on;
plot(range,real(savedtraj(ind+dimuhat*2,:)),'b','linewidth',1);hold on;
% xlim([min(range),max(range)])
xlim([1,400])
setgca(14);title(['(g) Mode (', num2str(ix),' , ', num2str(iy), ')'],'FontSize',13);
xlabel('t','fontsize',16)
lgnd=legend('Truth','Recovered with exact parameters','Recovered with estimated parameters');
set(lgnd,'FontSize',12);


subplot(4,3,1+3*3)
data1=savedtraj(1:dimuhat,timeplot2);
[ux,uy,xx,yy]=computeGBvel(rk,kk,L1,Dim_Grid,data1);

l = streamslice(xx,yy,ux,uy);
axis tight
set(l,'LineWidth',1);%xlabel('x','fontsize',16);%ylabel('y','fontsize',16);set(get(gca,'YLabel'),'Rotation',0)
set(l,'Color','k');
setgca(14);
ylabel(['(IV) t=',num2str(timeplot2)],'FontSize',13);set(get(gca,'YLabel'),'Rotation',0)
title('(j) Truth','fontsize',16)



subplot(4,3,3+3*2)
data1=savedtraj(1+2*dimuhat:3*dimuhat,timeplot1);
[ux,uy]=computeGBvel(rk,kk,L1,Dim_Grid,data1);
l = streamslice(xx,yy,ux,uy);
axis tight
set(l,'LineWidth',1);%xlabel('x','fontsize',16);%ylabel('y','fontsize',16);set(get(gca,'YLabel'),'Rotation',0)
set(l,'Color','b');
setgca(14);
title('(i) Estimated with L=60','fontsize',16)

subplot(4,3,3+3*3)
data1=savedtraj(1+2*dimuhat:3*dimuhat,timeplot2);
[ux,uy]=computeGBvel(rk,kk,L1,Dim_Grid,data1);
l = streamslice(xx,yy,ux,uy);
axis tight
set(l,'LineWidth',1);%xlabel('x','fontsize',16);%ylabel('y','fontsize',16);set(get(gca,'YLabel'),'Rotation',0)
set(l,'Color','b');
setgca(14);
title('(l) Estimated with L=60','fontsize',16)
