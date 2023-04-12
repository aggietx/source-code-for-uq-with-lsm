close all;clear;
% rmstime;rmsall;ccall
L=36;p=1;
if p==1
    if L==12
load data_p1_L12_para1d.mat;rmst(1,:)=rmstime;rms(1,:)=rmsall;cc(1,:)=ccall;
load data_p1_L12_para2d.mat;rmst(4,:)=rmstime;rms(4,:)=rmsall;cc(4,:)=ccall;
load data_p1_L12_5D;rmst(2,:)=rmstime;rms(2,:)=rmsall;cc(2,:)=ccall;
load data_p1_exact5D;rmst(5,:)=rmstime;rms(5,:)=rmsall;cc(5,:)=ccall;
load data_p1_L12_enkbf;rmst(3,:)=rmstime;rms(3,:)=rmsall;cc(3,:)=ccall;
    else
   load data_p1_L36_para1d.mat;rmst(1,:)=rmstime;rms(1,:)=rmsall;cc(1,:)=ccall;
load data_p1_L36_para2d.mat;rmst(4,:)=rmstime;rms(4,:)=rmsall;cc(4,:)=ccall;
load data_p1_L36_5D;rmst(2,:)=rmstime;rms(2,:)=rmsall;cc(2,:)=ccall;
load data_p1_exact5D;rmst(5,:)=rmstime;rms(5,:)=rmsall;cc(5,:)=ccall;
load data_p1_L36_enkbf;rmst(3,:)=rmstime;rms(3,:)=rmsall;cc(3,:)=ccall;
 
    end
else
        if L==12
load data_p05_L12_para1d.mat;rmst(1,:)=rmstime;rms(1,:)=rmsall;cc(1,:)=ccall;
load data_p05_L12_para2d.mat;rmst(4,:)=rmstime;rms(4,:)=rmsall;cc(4,:)=ccall;
load data_p05_L12_5D;rmst(2,:)=rmstime;rms(2,:)=rmsall;cc(2,:)=ccall;
load data_p05_exact5D;rmst(5,:)=rmstime;rms(5,:)=rmsall;cc(5,:)=ccall;
load data_p05_L12_enkbf;rmst(3,:)=rmstime;rms(3,:)=rmsall;cc(3,:)=ccall;
    else
   load data_p05_L36_para1d.mat;rmst(1,:)=rmstime;rms(1,:)=rmsall;cc(1,:)=ccall;
load data_p05_L36_para2d.mat;rmst(4,:)=rmstime;rms(4,:)=rmsall;cc(4,:)=ccall;
load data_p05_L36_5D;rmst(2,:)=rmstime;rms(2,:)=rmsall;cc(2,:)=ccall;
load data_p05_exact5D;rmst(5,:)=rmstime;rms(5,:)=rmsall;cc(5,:)=ccall;
load data_p05_L36_enkbf;rmst(3,:)=rmstime;rms(3,:)=rmsall;cc(3,:)=ccall;
 
         end
end
T=length(rmst);
np=length(cc);
figure()
subplot(3,1,1)
plot(1:T,rmst(1,:),'g','linewidth',1.5);hold on
plot(1:T,rmst(4,:),'b','linewidth',1.5);hold on
plot(1:T,rmst(2,:),'r','linewidth',1.5);hold on
plot(1:T,rmst(3,:),'k','linewidth',1.5);
plot(1:T,rmst(5,:),'c','linewidth',1.5);
% lgnd=legend('1D model','2D model','EnKBF Reduced','EnKBF perfect','Exact Reduced');
% set(lgnd,'fontsize',16;
setgca(14)
xlabel('t','fontsize',14)
xlim([1,T])
title('Mean RMS of velocity over fixed point','fontsize',14)


subplot(3,1,2)
plot(1:np,rms(1,:),'g','linewidth',1.5);hold on
plot(1:np,rms(4,:),'b','linewidth',1.5);hold on
plot(1:np,rms(2,:),'r','linewidth',1.5);hold on
plot(1:np,rms(3,:),'k','linewidth',1.5);
plot(1:np,rms(5,:),'c','linewidth',1.5);
% lgnd=legend('1D model','2D model','EnKBF Reduced','EnKBF perfect','Exact Reduced');
% set(lgnd,'fontsize',16;
setgca(14)
% xlabel('t','fontsize',14)
xlim([1,np])
title('Mean RMS of velocity over fixed point','fontsize',14)


subplot(3,1,3)
plot(1:np,cc(1,:),'g','linewidth',1.5);hold on
plot(1:np,cc(4,:),'b','linewidth',1.5);hold on
plot(1:np,cc(2,:),'r','linewidth',1.5);hold on
plot(1:np,cc(3,:),'k','linewidth',1.5);
plot(1:np,cc(5,:),'c','linewidth',1.5);
lgnd=legend('1D model','2D model','EnKBF Reduced','EnKBF perfect','Exact Reduced');
set(lgnd,'FontSize',14);
setgca(14)
xlabel('index of grid point','fontsize',14)
xlim([1,np])
title('Mean CC of velocity over fixed point','fontsize',14)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save uhate uhate;save estuhat estuhat;
% get  uhate.mat estuhat.mat

