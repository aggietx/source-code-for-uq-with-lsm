clear;
close all;
range=1:400;
plotcase=1;
load data_GB_rand3_3;ix=0;iy=1;ind=getind(ix,iy,kk);
figure();
subplot(3,1,1);
plot(range,savedtrajgb(ind,:),'k','linewidth',1.3);hold on;
if plotcase==1
plot(range,savedtrajgb(ind+dimuhat,:),'g','linewidth',1.3);hold on;
load data_GB_rand3_48;ix=0;iy=1;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'r','linewidth',1.3);hold on;
else
load data_GB_rand6_6;ix=0;iy=1;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'g','linewidth',1.3);hold on;
load data_GB_rand6_48;ix=0;iy=1;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'r','linewidth',1.3);hold on;
end
load data_GB_rand48_48;ix=0;iy=1;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'b','linewidth',1.3);hold on;
setgca(16);
xlim([1,400]);

subplot(3,1,2);
load data_GB_rand3_3;ix=2;iy=2;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind,:),'k','linewidth',1.3);hold on;
if plotcase==1
plot(range,savedtrajgb(ind+dimuhat,:),'g','linewidth',1.3);hold on;
load data_GB_rand3_48;ix=2;iy=2;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'r','linewidth',1.3);hold on;
else
load data_GB_rand6_6;ix=2;iy=2;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'g','linewidth',1.3);hold on;
load data_GB_rand6_48;ix=2;iy=2;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'r','linewidth',1.3);hold on;
end
load data_GB_rand48_48;ix=2;iy=2;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'b','linewidth',1.3);hold on;
xlim([1,400]);
setgca(16);

subplot(3,1,3);
load data_GB_rand3_3;ix=3;iy=3;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind,:),'k','linewidth',1.3);hold on;
if plotcase==1
plot(range,savedtrajgb(ind+dimuhat,:),'g','linewidth',1.3);hold on;
load data_GB_rand3_48;ix=3;iy=3;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'r','linewidth',1.3);hold on;
else
load data_GB_rand6_6;ix=3;iy=3;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'g','linewidth',1.3);hold on;
load data_GB_rand6_48;ix=3;iy=3;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'r','linewidth',1.3);hold on;
end
load data_GB_rand48_48;ix=3;iy=3;ind=getind(ix,iy,kk);
plot(range,savedtrajgb(ind+dimuhat,:),'b','linewidth',1.3);hold on;
setgca(16);

xlabel('t','fontsize',16)
xlim([1,400]);
if plotcase==1

lgnd=legend('Truth','3/3','3/48','48/48');
set(lgnd,'FontSize',14);
else
    lgnd=legend('Truth','6/6','6/48','48/48');
end