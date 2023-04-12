%% plot pdf
clear;
figure();%% all are rg data
% load trajcase1;
% load trajcase2;
% load traj2;%%%% p=2;
% load traj1;%%%% p=2;
% load traj;%%%% p=2;
load trajc1;traj=real(traj(1:3,:));
Naut=500;
Dt=0.01;
for ii=1:3
data=traj(ii,:);

subplot(3,3,3*ii-2);
plot(Dt:Dt:Dt*length(data),data,'b','linewidth',1);
xlim([Dt,Dt*length(data)])
setgca(14);
% if ii==1
%     ylabel('u','Fontsize',12)
% elseif ii==2
%      ylabel('\psi_1','Fontsize',12)
% else
%     ylabel('\psi_2','Fontsize',12)
% end

if ii==1
    ylabel('(0,0)','Fontsize',12)
elseif ii==2
     ylabel('(1,1)','Fontsize',12)
else
    ylabel('(2,2)','Fontsize',12)
end

set(get(gca,'YLabel'),'Rotation',0)
np=400;
data=real(data);
data=data(:);
mu=mean(data);sig=var(data);
[~,r1]=ksdensity(data);
maxdata=max(r1);mindata=min(r1);
dx=(maxdata-mindata)/np;
xx=mindata:dx:maxdata;
[prob,r1]=ksdensity(data,xx);
%% gaussfit
datan=mu+sqrt(sig)*randn(length(data),1);
% [probn,rn]=ksdensity(datan);
maxdatan=max(datan);mindatan=min(datan);
dx=(maxdatan-mindatan)/np;
xxn=mindatan:dx:maxdatan;
  [probn,rn]=ksdensity(datan,xxn);
subplot(3,3,3*ii-1);
   plot(r1,prob,'blue','Linewidth',2);
hold on;plot(rn,probn,'--r','Linewidth',2);
xlim([mindatan, maxdatan]);setgca(14)
if ii==3
  lgnd=legend('PDF',...
    'Gaussian fit');
set(lgnd,'FontSize',12);
end
subplot(3,3,3*ii);
%    plot(r1,log(prob),'blue','Linewidth',2);
% hold on;plot(rn,log(probn),'--r','Linewidth',2);
% xlim([mindatan, maxdatan]);

acf1 = autocorr(data,Naut);
plot(0:Dt:Naut*Dt,acf1,'r','linewidth',2);
xlim([0,Naut*Dt]);



% title('PDF','Fontsize',16)

setgca(14)


% title('(II) PDF','Fontsize',16)
end