close all;
L=12;p=1;
if p==1
    if L==12
load data_p1_L12_para1d.mat;data(1,:)=uhate(1,:);data1(1,:)=uhate(2,:);data2(1,:)=uhate(4,:);
data(2,:)=estuhat(1,:);data1(2,:)=estuhat(2,:);data2(2,:)=estuhat(4,:);
load data_p1_L12_para2d.mat;data(5,:)=estuhat(1,:);
load data_p1_L12_5D;data(3,:)=uhate(1,:);data1(3,:)=uhate(2,:);data2(3,:)=uhate(4,:);
load data_p1_L12_enkbf;data(4,:)=uhate(1,:);data1(4,:)=uhate(2,:);data2(4,:)=uhate(4,:);
    else
   load data_p1_L36_para1d.mat;data(1,:)=uhate(1,:);data1(1,:)=uhate(2,:);data2(1,:)=uhate(4,:);
   data(2,:)=estuhat(1,:);data1(2,:)=estuhat(2,:);data2(2,:)=estuhat(4,:);
load data_p1_L36_para2d.mat;data(5,:)=estuhat(1,:);
load data_p1_L36_5D;data(3,:)=uhate(1,:);data1(3,:)=uhate(2,:);data2(3,:)=uhate(4,:);
load data_p1_L36_enkbf;data(4,:)=uhate(1,:);data1(4,:)=uhate(2,:);data2(4,:)=uhate(4,:);
 
    end
else
        if L==12
load data_p05_L12_para1d.mat;data(1,:)=uhate(1,:);data1(1,:)=uhate(2,:);data2(1,:)=uhate(4,:);
data(2,:)=estuhat(1,:);data1(2,:)=estuhat(2,:);data2(2,:)=estuhat(4,:);
load data_p05_L12_para2d.mat;data(5,:)=estuhat(1,:);
load data_p05_L12_5D;data(3,:)=uhate(1,:);data1(3,:)=uhate(2,:);data2(3,:)=uhate(4,:);
load data_p05_L12_enkbf;data(4,:)=uhate(1,:);data1(4,:)=uhate(2,:);data2(4,:)=uhate(4,:);
    else
   load data_p05_L36_para1d.mat;data(1,:)=uhate(1,:);data1(1,:)=uhate(2,:);data2(1,:)=uhate(4,:);
   data(2,:)=estuhat(1,:);data1(2,:)=estuhat(2,:);data2(2,:)=estuhat(4,:);
load data_p05_L36_para2d.mat;data(5,:)=estuhat(1,:);
load data_p05_L36_5D;data(3,:)=uhate(1,:);data1(3,:)=uhate(2,:);data2(3,:)=uhate(4,:);
load data_p05_L36_enkbf;data(4,:)=uhate(1,:);data1(4,:)=uhate(2,:);data2(4,:)=uhate(4,:);
 
         end
end
data=real(data);
data1=real(data1);data2=real(data2);
 T=size(data,2);
figure()
subplot(3,1,1)
plot(1:T,data(1,:),'r','linewidth',1.5);hold on
plot(1:T,real(data(2,:)),'b','linewidth',1.5);hold on
plot(1:T,real(data(5,:)),'c','linewidth',1.5);hold on
plot(1:T,data(3,:),'g','linewidth',1.5);hold on
plot(1:T,real(data(4,:)),'k','linewidth',1.5);

lgnd=legend('Exact','1D model','2D model','EnKBF Reduced','EnKBF perfect');
set(lgnd,'FontSize',10);
setgca(16)
% xlabel('t','fontsize',24)
xlim([1,T])
title('u','fontsize',14)

subplot(3,1,2)
plot(1:T,data1(1,:),'r','linewidth',1.5);hold on
plot(1:T,real(data1(2,:)),'b','linewidth',1.5);hold on
plot(1:T,data1(3,:),'g','linewidth',1.5);hold on
plot(1:T,real(data1(4,:)),'k','linewidth',1.5);hold on
% plot(1:T,real(data(5,:)),'c','linewidth',1.5);hold on
% lgnd=legend('Exact','1D model','EnKBF 5D','EnKBF perfect','2D');
% set(lgnd,'FontSize',12);
setgca(16)
% xlabel('t','fontsize',24)
xlim([1,T])
title('\psi_1','fontsize',14)

subplot(3,1,3)
plot(1:T,data2(1,:),'r','linewidth',1.5);hold on
plot(1:T,real(data2(2,:)),'b','linewidth',1.5);hold on
plot(1:T,data2(3,:),'g','linewidth',1.5);hold on
plot(1:T,real(data2(4,:)),'k','linewidth',1.5);hold on
% plot(1:T,real(data(5,:)),'c','linewidth',1.5);hold on
% lgnd=legend('Exact','1D model','EnKBF 5D','EnKBF perfect','2D');
% set(lgnd,'FontSize',12);
setgca(16)
xlabel('t','fontsize',24)
xlim([1,T])
title('\psi_2','fontsize',14)