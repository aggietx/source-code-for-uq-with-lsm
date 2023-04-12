clear 
% close all;
%%%% clearlay;save  data2d_layer1_L12_p2
% % clearlay;save data1d_layer1_L3_p2
% % clearlay;save data1dreduced_layer1_L3_p2
% % clearlay;save data_enkbf1_L3_p2
% % clearlay;save data_enkbf1_reduced_L3_p2
plotvel=0;S=1.5;
% % load data2d_layer1_L12
% % load data1dreduced_layer1_L12
% % load data_enkbf1_L12
% % load data_enkbf1_reduced_L12 

% load data1d_layer1_L12
load data1d_layer1_L6_p2
time1=find(uexact==max(uexact));
time2=find(abs(uexact-max(uexact)/2)<0.2);
time2=time2(6);
time3=find(abs(uexact)<0.2);
time3=time3(6);


figure();

subplot(3,3,1);
[velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,uhate(:,time1));
if plotvel
    quiver(xx,yy,velx,vely,S, 'linewidth',1.5,'color','k')
else

l = streamslice(xx,yy,velx,vely);
axis tight
set(l,'LineWidth',1.5)
set(l,'Color','k');

end
xlim([0,2*pi]);ylim([0,2*pi]);setgca(14);
ylabel('(I) Truth','FontSize',16);set(get(gca,'Ylabel'),'Rotation',0)
title(['(a) t=',num2str(time1)],'FontSize',16);

subplot(3,3,4);
% load data_enkbf1_L12
load data_enkbf1_L6_p2
[velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,estuhat(:,time1));
if plotvel
    quiver(xx,yy,velx,vely,S, 'linewidth',1.5,'color','r')
else

l = streamslice(xx,yy,velx,vely);
axis tight
set(l,'LineWidth',1.5)
set(l,'Color','r');

end
xlim([0,2*pi]);ylim([0,2*pi]);setgca(14);
ylabel('(II) EnKBF all modes','FontSize',16);set(get(gca,'Ylabel'),'Rotation',0)

subplot(3,3,7);
% load data1d_layer1_L12
load data1d_layer1_L6_p2
[velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,estuhat(:,time1));
if plotvel
    quiver(xx,yy,velx,vely,S, 'linewidth',1.5,'color','b')
else

l = streamslice(xx,yy,velx,vely);
axis tight
set(l,'LineWidth',1.5)
set(l,'Color','b');

end
xlim([0,2*pi]);ylim([0,2*pi]);setgca(14);
ylabel('(III) 1D Estimation','FontSize',16);set(get(gca,'Ylabel'),'Rotation',0)
% xlabel('x','FontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,2);
[velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,uhate(:,time2));
if plotvel
    quiver(xx,yy,velx,vely,S, 'linewidth',1.5,'color','k')
else

l = streamslice(xx,yy,velx,vely);
axis tight
set(l,'LineWidth',1.5)
set(l,'Color','k');

end
xlim([0,2*pi]);ylim([0,2*pi]);setgca(14);
% ylabel('(I) Truth','FontSize',16);set(get(gca,'Ylabel'),'Rotation',0)
title(['(b) t=',num2str(time2)],'FontSize',16);

subplot(3,3,5);
% load data_enkbf1_L12
load data1d_layer1_L6_p2
[velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,estuhat(:,time2));
if plotvel
    quiver(xx,yy,velx,vely,S, 'linewidth',1.5,'color','r')
else

l = streamslice(xx,yy,velx,vely);
axis tight
set(l,'LineWidth',1.5)
set(l,'Color','r');

end
xlim([0,2*pi]);ylim([0,2*pi]);setgca(14);
% ylabel('(II) EnKBF','FontSize',16);set(get(gca,'Ylabel'),'Rotation',0)

subplot(3,3,8);
% load data1d_layer1_L12
load data1d_layer1_L6_p2
[velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,estuhat(:,time2));
if plotvel
    quiver(xx,yy,velx,vely,S, 'linewidth',1.5,'color','b')
else

l = streamslice(xx,yy,velx,vely);
axis tight
set(l,'LineWidth',1.5)
set(l,'Color','b');

end
xlim([0,2*pi]);ylim([0,2*pi]);setgca(14);
% ylabel('(III) 1D Estimation','FontSize',16);set(get(gca,'Ylabel'),'Rotation',0)
% xlabel('x','FontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,3);
[velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,uhate(:,time3));
if plotvel
    quiver(xx,yy,velx,vely,S, 'linewidth',1.5,'color','k')
else

l = streamslice(xx,yy,velx,vely);
axis tight
set(l,'LineWidth',1.5)
set(l,'Color','k');

end
xlim([0,2*pi]);ylim([0,2*pi]);setgca(14);
% ylabel('(I) Truth','FontSize',16);set(get(gca,'Ylabel'),'Rotation',0)
title(['(c) t=',num2str(time3)],'FontSize',16);

subplot(3,3,6);
% load data_enkbf1_L12
load data_enkbf1_L6_p2
[velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,estuhat(:,time3));
if plotvel
    quiver(xx,yy,velx,vely,S, 'linewidth',1.5,'color','r')
else

l = streamslice(xx,yy,velx,vely);
axis tight
set(l,'LineWidth',1.5)
set(l,'Color','r');

end
xlim([0,2*pi]);ylim([0,2*pi]);setgca(14);
% ylabel('(II) EnKBF','FontSize',16);set(get(gca,'Ylabel'),'Rotation',0)

subplot(3,3,9);
% load data1d_layer1_L12
load data1d_layer1_L6_p2
[velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,estuhat(:,time3));
if plotvel
    quiver(xx,yy,velx,vely,S, 'linewidth',1.5,'color','b')
else

l = streamslice(xx,yy,velx,vely);
axis tight
set(l,'LineWidth',1.5)
set(l,'Color','b');

end
xlim([0,2*pi]);ylim([0,2*pi]);setgca(14);
% ylabel('(III) 1D Estimation','FontSize',16);set(get(gca,'Ylabel'),'Rotation',0)
% xlabel('x','FontSize',16);