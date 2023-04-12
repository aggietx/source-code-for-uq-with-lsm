clear;close all
rng(10);
load u_hat;T=400;u_hat=u_hat*10^5;
dt=0.005;Kmax=5;K_max=5;
[ky,kx]=meshgrid([0:Kmax,-Kmax:-1],[0:Kmax,-Kmax:-1]);
kk=[kx(:),ky(:)]';

funfit=0;maxlag=500;Naut=maxlag;%% for plot
datafit=u_hat;
% datafit=gamma_mean_trace;
% datafit=gamma_mean_trace_smoother;
% datafit=gamma_mean_trace_smoother_sample;
[dk,omegak,mu,sigmak2]=fit_ou_fun(funfit,K_max,maxlag,dt,datafit,kk);
T1=T*5;N1=T/dt;
uhatfit=zeros(size(u_hat,1),N1);
uhatfit1=zeros(1,N1);
sigmafit=form_sigmatrix1((2*Kmax+1)^2,sigmak2,kk,Kmax);
disp('free run');
for i=2:N1
uhatfit(:,i)=uhatfit(:,i-1)+(-dk+1i*omegak).*uhatfit(:,i-1)*dt+mu*dt+sigmafit*randn(size(u_hat,1),1)*sqrt(dt);
% uhatfit1(:,i)=uhatfit1(:,i-1)+(-dk(25)+1i*omegak(25)).*uhatfit1(:,i-1)*dt+mu(25)*dt+...
%     sigmak2(25)*(randn*sqrt(dt)+1i*randn*sqrt(dt));

end
ix=2;iy=2;
[a]=find(kk(1,:)==ix);
[b]=find(kk(2,:)==iy);
ind=intersect(a,b);
data2=real(u_hat(ind,1:1:end));%%%% exact data
data1=real(uhatfit(ind,1:end));
 data1=data1(:);data2=data2(:);

 [prob,r1]=ksdensity(data1);
figure();
% subplot(2,1,1)
plot(r1,prob,'r','Linewidth',2); set(gca,'fontsize',16)
hold on
 [prob,r1]=ksdensity(data2);
plot(r1,prob,'g','Linewidth',2); set(gca,'fontsize',16)
lgnd=legend('truth', 'gauss fit');
set(lgnd,'FontSize',12,'Location', 'Best');
title([' real, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
% return
% filenameStrPdf = sprintf('pdfou_%d_%d.pdf',ix,iy); 
% print(gcf, '-dpdf', '-r600', filenameStrPdf)
% filenameStrfig = sprintf('pdfou_%d_%d.fig',ix,iy); 
% savefig(filenameStrfig)




%%%


acf1 = autocorr(data1,Naut);
acf2 = autocorr(data2,Naut);

figure()
plot(0:dt:Naut*dt,acf1,'r','linewidth',2);hold on;
plot(0:dt:Naut*dt,acf2,'g','linewidth',2);hold on;
set(gca,'fontsize',16)
lgnd=legend('fit free run', ' truth');
set(lgnd,'FontSize',14,'Location', 'Best');
xlabel('t','fontsize',16);

return
title([' real u, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
filenameStrPdf = sprintf('acfu_%d_%d.pdf',ix,iy);   
print(gcf, '-dpdf', '-r600',filenameStrPdf)
filenameStrfig = sprintf('acfu_%d_%d.fig',ix,iy); 
savefig(filenameStrfig)
