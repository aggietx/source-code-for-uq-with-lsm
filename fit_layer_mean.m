%%%% only test mean field
%%%% run fit_layered first to generate the fit data
maxlag0=100000;funfit=0;ts=1;fprintf('Lag are  %d \n',maxlag0);
fprintf('testcase is %d\n',testcase);
fprintf('funfit is %d\n',funfit);
dkfit0=zeros(1,1);omegakfit0=zeros(1,1);
Ffit0=zeros(1,1);sigmak20=zeros(1,1);
[dkfit0(1),omegakfit0(1),Ffit0(1),sigmak20(1)]=fit_ou_fun(funfit,Kmax,maxlag0,dt,traj1(1,ts/dt:end),kk1dfullfake(:,1));
sigmafit0=sigmak20;


disp('free run..')
dim=length(dkfit0);
ufree=zeros(dim, N);
uold=zeros(dim,1);
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end    
   u1=uold+(-dkfit0+1i*omegakfit0).*uold*dt+Ffit0*dt+sqrt(dt)*sigmafit0*randn(dim,1);
   uold=u1;
   ufree(:,i)=u1;
end
ufree1=ufree;
clear ufree
% id=2;
% [prob,r1]=plotpdf(real(traj1(id,ts/dt:end)));b=[prob;r1];save b b
% % plotpdf(real(ufree(id,ts/dt:end)));
% data=real(ufree(id,ts/dt:end));[proba,r1a]=plotpdf(data);a=[proba;r1a];save a a
% 
% load a;load b;close all
% figure();plot(b(2,:),b(1,:),'k');hold on
% plot(a(2,:),a(1,:),'b')

range=10/dt:T/dt;Naut=5000;aa=100;

for i=1:1
    figure()
    subplot(2,1,1)
    [prob,r1]=ksdensity(real(traj(i,range)));
    plot(r1,prob,'k','Linewidth',1.5);
    [prob,r1]=ksdensity(real(ufree1(i,range)));hold on
    plot(r1,prob,'b','Linewidth',1.5);
    setgca(16)
    
   subplot(2,1,2)
    if i==1
    acf1=autocorr(real(traj(i,range)),Naut*aa);
    plot(0:dt:Naut*aa*dt,acf1,'k','Linewidth',1.5);
    acf1=autocorr(real(ufree1(i,range)),Naut*aa);hold on
    plot(0:dt:Naut*aa*dt,acf1,'b','Linewidth',1.5);    
    else
     acf1=autocorr(real(traj(i,range)),Naut);
    plot(0:dt:Naut*dt,acf1,'k','Linewidth',1.5); 
    acf1=autocorr(real(ufree1(i,range)),Naut);hold on
    plot(0:dt:Naut*dt,acf1,'b','Linewidth',1.5);    
    end    
    setgca(16)

  lgnd=legend('Truth',...
    'Free run fit');
set(lgnd,'FontSize',12);
 filenameStrPdf = sprintf('oufitcase%dmode%d.pdf',testcase,i-1);   
print(gcf, '-dpdf', '-r600', filenameStrPdf)
end


