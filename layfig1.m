% close all;
clear
% load  trajreg1;p=1;
% load  trajreg2;p=1;
load traj1;trajreg1=traj(:,1:20000);
ene1=getene(traj);ene1=ene1(2:end);
load traj2;trajreg2=traj(:,1:20000);
ene2=getene(traj);ene2=ene2(2:end);

% ene1=[7.3964
%    27.6042
%     3.8815
%     2.8522
%     2.2287
%     1.6376];
% ene2=[1.2393
%     5.1077
%     0.9189
%     0.6539
%     0.6332
%     0.4771];
Dt=0.1;N=length(trajreg1);T=Dt*N;
Tr=Dt:Dt:T;Kmax=size(trajreg1,1)-1;
nx=20;x=0:2*pi/nx:2*pi;
uy1d1=zeros(length(x),N/10);
uy1d2=zeros(length(x),N/10);
for j=1:N/10
uy1d1(:,j)=computeuy1d(trajreg1(2:end,j*10),Kmax,x');
uy1d2(:,j)=computeuy1d(trajreg2(2:end,j*10),Kmax,x');
end

figure();

subplot(4,8,[1,2,3])
plot(Tr,trajreg1(1,:),'b','linewidth',1)
setgca(13);
title('(a) Trajectory (Regime I)','FontSize',14)
% %ylabel('(I) u','FontSize',14);set(get(gca,'%ylabel'),'Rotation',0)

subplot(4,8,4)
setgca(13);title('(b) PDF','FontSize',14)

subplot(4,8,[5,6,7])
plot(Tr,trajreg2(1,:),'b','linewidth',1);setgca(13);
title('(c) Trajectory (Regime II)','FontSize',14)
% %ylabel('u','FontSize',14);set(get(gca,'%ylabel'),'Rotation',0)

subplot(4,8,8);
title('(d) PDF','FontSize',14)
setgca(13);
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,8,[1+8,2+8,3+8])
plot(Tr,real(trajreg1(2,:)),'b','linewidth',1)
setgca(13);%ylabel('(II) \psi_1','FontSize',14);set(get(gca,'%ylabel'),'Rotation',0)
% xlabel('t','FontSize',14);

subplot(4,8,4+8)
setgca(13);

subplot(4,8,[5+8,6+8,7+8])
plot(Tr,real(trajreg2(2,:)),'b','linewidth',1)
setgca(13);

subplot(4,8,8+8)
setgca(13);
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,8,[1+8*2,2+8*2,3+8*2])
plot(Tr,real(trajreg1(3,:)),'b','linewidth',1)
setgca(13);%ylabel('(III) \psi_2','FontSize',14);set(get(gca,'%ylabel'),'Rotation',0)

subplot(4,8,4+8*2)
setgca(13);

subplot(4,8,[5+8*2,6+8*2,7+8*2]);setgca(13);
plot(Tr,real(trajreg2(3,:)),'b','linewidth',1)

subplot(4,8,8+8*2)
setgca(13);
lgnd=legend('Truth','Gauss fit');
set(lgnd,'FontSize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,8,[1+8*3,2+8*3,3+8*3])
% plot(Tr,real(trajreg1(2,:)),'b','linewidth',1)
imagesc([1, T],[0,2*pi],uy1d1);colorbar
xlabel('t','FontSize',14);
colorbar
% %ylabel('(IV) x','FontSize',14);set(get(gca,'%ylabel'),'Rotation',0)
set(gca,'ydir','normal')
setgca(13);

subplot(4,8,4+8*2)
setgca(13);

subplot(4,8,[5+8*3,6+8*3,7+8*3])
plot(Tr,real(trajreg2(2,:)),'b','linewidth',1)
imagesc([1, T],[0,2*pi],uy1d2);
set(gca,'ydir','normal')
colorbar
% colormap(cmocean('curl'))
xlabel('t','FontSize',14);

setgca(13);
% subplot(4,8,8+8*2)
%%%%%%
subplot(4,8,4+8*3)
plot(1:Kmax,ene1,'*-b','linewidth',2);setgca(13);
xlabel('k','FontSize',14);
xlim([1,Kmax]);
title('(e) Energy','FontSize',13);
%%%%%%
subplot(4,8,8+8*3)
plot(1:Kmax,ene2,'*-b','linewidth',2);setgca(13);
xlabel('k','FontSize',14);
xlim([1,Kmax]);
title('(f) Energy','FontSize',13);
%%%%% 
