% close all;
pid=5;full=1;puregb=1;
if 0
[dke,omegake,mue,sigmak2e]=fit_ou_fun_sw(funfit,K_max,maxlag1,maxlag2,dt,u_hat(:,tfit/dt:end),kk);
alldkfit=[alldkfit,dke];allsigmak2=[allsigmak2,sigmak2e];
save alldkfit alldkfit;save allomegakfit allomegakfit;save allsigmak2 allsigmak2
get alldkfit.mat  allsigmak2.mat allomegakfit.mat
end
load alldkfit;load allsigmak2;load allomegakfit
if full
K_max=5;redind=(2:(2*K_max+1)^2)';
% K_max=3;redind=(1+(2*K_max+1)^2:3*(2*K_max+1)^2)';
% K_max=3;redind=(2:(2*K_max+1)^2)';
else
redind=[];
[ky,kx]=meshgrid([0:K_max,-K_max:-1],[0:K_max,-K_max:-1]);
kk=[kx(:),ky(:)]';
for ii=1:K_max
  redind=[redind;getind(ii,ii,kk)];  
end
for ii=1:K_max
  redind=[redind;getind(ii,0,kk)];  
end
end

%% plot estimation comparison
if puregb
    figure();subplot(2,1,1);
else
figure();subplot(3,1,1);
end
% plot(-dkexact(redind),'*-r','linewidth',2);hold on;
plot(-alldkfit(redind,end-1),'*-r','linewidth',2);hold on;
plot(alldkfit(redind,1),'*-b','linewidth',2);hold on;
plot(alldkfit(redind,pid+1),'*-g','linewidth',2);
plot(alldkfit(redind,end),'*-k','linewidth',2);
% lgnd=legend('Exact','Initial','Estimated');
% set(lgnd,'FontSize',18);
setgca(18)
title('d_k','FontSize',18);
xlim([1,length(redind)])

if puregb
    subplot(2,1,2);
else
subplot(3,1,3);
end
% plot(real(sigmak2exact(redind))/sqrt(2),'*-r','linewidth',2);hold on;
plot(real(allsigmak2((redind),end-1)),'*-r','linewidth',2);hold on
plot(real(allsigmak2((redind),1)),'*-b','linewidth',2);hold on
plot(real(allsigmak2((redind),pid+1)),'*-g','linewidth',2);
plot(real(allsigmak2((redind),end)),'*-k','linewidth',2);hold on
% lgnd=legend('Exact','Initial','Estimated');
% lgnd=legend('Exact','Initial','Estimated','Exact uhat fit');
% set(lgnd,'FontSize',18);
if puregb
    lgnd=legend('Exact','Initial','Estimated','Exact uhat fit');
set(lgnd,'FontSize',15);
end
setgca(18)
title('\sigma_k','FontSize',18);
xlim([1,length(redind)])


% figure();
if puregb~=1
subplot(3,1,2);
plot(real(allomegakfit((redind),end-1)),'*-r','linewidth',2);hold on
plot(real(allomegakfit((redind),1)),'*-b','linewidth',2);hold on
plot(real(allomegakfit((redind),pid+1)),'*-g','linewidth',2);
plot(real(allomegakfit((redind),end)),'*-k','linewidth',2);hold on
% lgnd=legend('Exact','Initial','Estimated');
lgnd=legend('Exact','Initial','Estimated','Exact uhat fit');
set(lgnd,'FontSize',15);
setgca(18)
title('\omega_k','FontSize',18);
xlim([1,length(redind)])
end