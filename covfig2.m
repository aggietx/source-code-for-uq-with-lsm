clear;
close all;

figure()
%%%%%%

load datagbcase1_L12;range=1:T;ix=1;iy=0;ind=getind(ix,iy,kk);data=real(savedtrajgb);
subplot(3,3,1)
plot(range,data(ind,:),'k','linewidth',1.5);hold on;
plot(range,data(ind+dimuhat,:),'b','linewidth',1.5);hold on;
load datagbcase2_L12;ix=1;iy=0;ind=getind(ix,iy,kk);data=real(savedtrajgb);
% plot(range,data(ind+dimuhat,:),'r','linewidth',1.5);hold on;
load datagbcase3_L12;ix=1;iy=0;
ind=getind(ix,iy,kk);data=real(savedtrajgb);
plot(range,data(ind+dimuhat,:),'g','linewidth',1.5);hold on;
title('(a) GB','FontSize',13);setgca(13)
ylabel('(I) (1,0),','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)

%%%%%%%

load datagbcase1_L12;ix=2;iy=0;ind=getind(ix,iy,kk);data=real(savedtrajgb);
subplot(3,3,4)
plot(range,data(ind,:),'k','linewidth',1.5);hold on;
plot(range,data(ind+dimuhat,:),'b','linewidth',1.5);hold on;
load datagbcase2_L12;ix=2;iy=0;ind=getind(ix,iy,kk);data=real(savedtrajgb);
% plot(range,data(ind+dimuhat,:),'r','linewidth',1.5);hold on;
load datagbcase3_L12;
ix=2;iy=0;ind=getind(ix,iy,kk);data=real(savedtrajgb);
plot(range,data(ind+dimuhat,:),'g','linewidth',1.5);hold on;
% title('(a) GB','FontSize',13);
setgca(13)
ylabel('(II) (2,0),','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)

%%%%%%%
load datagbcase1_L12;ix=2;iy=2;ind=getind(ix,iy,kk);data=real(savedtrajgb);
subplot(3,3,7)
plot(range,data(ind,:),'k','linewidth',1.5);hold on;
plot(range,data(ind+dimuhat,:),'b','linewidth',1.5);hold on;
load datagbcase2_L12;ix=2;iy=2;ind=getind(ix,iy,kk);data=real(savedtrajgb);
% plot(range,data(ind+dimuhat,:),'r','linewidth',1.5);hold on;
load datagbcase3_L12;
ix=2;iy=2;ind=getind(ix,iy,kk);data=real(savedtrajgb);
plot(range,data(ind+dimuhat,:),'g','linewidth',1.5);hold on;
% title('(a) GB','FontSize',13);
xlabel('t','FontSize',13);
setgca(13)
ylabel('(III) (2,2),','FontSize',13);set(get(gca,'YLabel'),'Rotation',0)

%%%%%%%%
load dataswcase1_L30;
time1=find(rmstimefilgr==min(rmstimefilgr));
% min(rmstimefilgr)
  time2=find(rmstimefilgr==max(rmstimefilgr)); 
%   max(rmstimefilgr)

range=1:T;ix=1;iy=0;ind=getind(ix,iy,kk);ind=ind(1);data=real(savedtrajsw);
subplot(3,3,2)
plot(range,data(ind,:),'k','linewidth',1.5);hold on;
plot(range,data(ind+dimuhat,:),'b','linewidth',1.5);hold on;
load dataswcase2_L30;ix=1;iy=0;ind=getind(ix,iy,kk);ind=ind(1);data=real(savedtrajsw);
% plot(range,data(ind+dimuhat,:),'r','linewidth',1.5);hold on;
load dataswcase3_L30;ix=1;iy=0;
ind=getind(ix,iy,kk);ind=ind(1);data=real(savedtrajsw);
plot(range,data(ind+dimuhat,:),'g','linewidth',1.5);hold on;

plot(time1,data(ind,time1),'sr','MarkerSize',5,'linewidth',2);hold on
plot(time2,data(ind,time2),'sr','MarkerSize',5,'linewidth',2);hold on
title('(b) GB of SW','FontSize',13);setgca(13)
%%%%%%
load dataswcase1_L30;
range=1:T;ix=1;iy=0;ind=getind(ix,iy,kk);ind=ind(2);data=real(savedtrajsw);
subplot(3,3,3)
plot(range,data(ind,:),'k','linewidth',1.5);hold on;
plot(range,data(ind+dimuhat,:),'b','linewidth',1.5);hold on;
load dataswcase2_L30;ix=1;iy=0;ind=getind(ix,iy,kk);ind=ind(2);data=real(savedtrajsw);
% plot(range,data(ind+dimuhat,:),'r','linewidth',1.5);hold on;
load dataswcase3_L30;ix=1;iy=0;
ind=getind(ix,iy,kk);ind=ind(2);data=real(savedtrajsw);
plot(range,data(ind+dimuhat,:),'g','linewidth',1.5);hold on;
plot(time1,data(ind,time1),'sr','MarkerSize',5,'linewidth',2);hold on
plot(time2,data(ind,time2),'sr','MarkerSize',5,'linewidth',2);hold on
title('(c) GR of SW','FontSize',13);setgca(13)
%%%%%%

%%%%%%%%
load dataswcase1_L30;
range=1:T;ix=2;iy=0;ind=getind(ix,iy,kk);ind=ind(1);data=real(savedtrajsw);
subplot(3,3,5)
plot(range,data(ind,:),'k','linewidth',1.5);hold on;
plot(range,data(ind+dimuhat,:),'b','linewidth',1.5);hold on;
load dataswcase2_L30;ix=2;iy=0;ind=getind(ix,iy,kk);ind=ind(1);data=real(savedtrajsw);
% plot(range,data(ind+dimuhat,:),'r','linewidth',1.5);hold on;
load dataswcase3_L30;ix=2;iy=0;
ind=getind(ix,iy,kk);ind=ind(1);data=real(savedtrajsw);
plot(range,data(ind+dimuhat,:),'g','linewidth',1.5);hold on;
% title('(b) GB of SW','FontSize',13);
setgca(13)
%%%%%%
load dataswcase1_L30;
range=1:T;ix=2;iy=0;ind=getind(ix,iy,kk);ind=ind(2);data=real(savedtrajsw);
subplot(3,3,6)
plot(range,data(ind,:),'k','linewidth',1.5);hold on;
plot(range,data(ind+dimuhat,:),'b','linewidth',1.5);hold on;
load dataswcase2_L30;ix=2;iy=0;ind=getind(ix,iy,kk);ind=ind(2);data=real(savedtrajsw);
% plot(range,data(ind+dimuhat,:),'r','linewidth',1.5);hold on;
load dataswcase3_L30;ix=2;iy=0;
ind=getind(ix,iy,kk);ind=ind(2);data=real(savedtrajsw);
plot(range,data(ind+dimuhat,:),'g','linewidth',1.5);hold on;
% title('(c) GR of SW','FontSize',13);
setgca(13)
%%%%%%
%%%%%%%%
load dataswcase1_L30;
range=1:T;ix=2;iy=2;ind=getind(ix,iy,kk);ind=ind(1);data=real(savedtrajsw);
subplot(3,3,8)
plot(range,data(ind,:),'k','linewidth',1.5);hold on;
plot(range,data(ind+dimuhat,:),'b','linewidth',1.5);hold on;
load dataswcase2_L30;ix=2;iy=2;ind=getind(ix,iy,kk);ind=ind(1);data=real(savedtrajsw);
% plot(range,data(ind+dimuhat,:),'r','linewidth',1.5);hold on;
load dataswcase3_L30;ix=2;iy=2;
ind=getind(ix,iy,kk);ind=ind(1);data=real(savedtrajsw);
plot(range,data(ind+dimuhat,:),'g','linewidth',1.5);hold on;
% title('(b) GB of SW','FontSize',13);
setgca(13)
xlabel('t','FontSize',13);
%%%%%%
load dataswcase1_L30;
range=1:T;ix=2;iy=2;ind=getind(ix,iy,kk);ind=ind(2);data=real(savedtrajsw);
subplot(3,3,9)
plot(range,data(ind,:),'k','linewidth',1.5);hold on;
plot(range,data(ind+dimuhat,:),'b','linewidth',1.5);hold on;
load dataswcase2_L30;ix=2;iy=2;ind=getind(ix,iy,kk);ind=ind(2);data=real(savedtrajsw);
% plot(range,data(ind+dimuhat,:),'r','linewidth',1.5);hold on;
load dataswcase3_L30;ix=2;iy=2;
ind=getind(ix,iy,kk);ind=ind(2);data=real(savedtrajsw);
plot(range,data(ind+dimuhat,:),'g','linewidth',1.5);hold on;
% title('(c) GR of SW','FontSize',13);
xlabel('t','FontSize',13);
setgca(13)
% lgnd=legend('Truth','Full R','Diagonal R','Constant R');
lgnd=legend('Truth','Full R','Constant R');
set(lgnd,'FontSize',12);

%%%%%%





