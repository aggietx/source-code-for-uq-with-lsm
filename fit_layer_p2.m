%%%% run fit_layered first to generate the fit data
 maxlag=4000;maxlag0=200000;funfit=1;ts=1;fprintf('Lag are %d %d \n',maxlag0,maxlag);
fprintf('testcase is %d\n',testcase);
fprintf('funfit is %d\n',funfit);
dkfit=zeros(2*Kmax+1,1);omegakfit=zeros(2*Kmax+1,1);
Ffit=zeros(2*Kmax+1,1);sigmak2=zeros(2*Kmax+1,1);
[dkfit(1),omegakfit(1),Ffit(1),sigmak2(1)]=fit_ou_fun(funfit,Kmax,maxlag0,dt,traj1(1,ts/dt:end),kk1dfullfake(:,1));
[dkfit(2:end),omegakfit(2:end),Ffit(2:end),sigmak2(2:end)]=fit_ou_fun(funfit,Kmax,maxlag,dt,traj1(2:end,ts/dt:end),kk1dfullfake(:,2:end));
sigmafit=form_sigmatrix1d(2*Kmax+1,sigmak2,kk1dfullfake,Kmax);
% return
if testcase==1
   save dkfitc1 dkfit
   save  omegakfitc1 omegakfit
   save  Ffitc1 Ffit
   save sigmafitc1 sigmafit
else
       save dkfitc2 dkfit
   save  omegakfitc2 omegakfit
   save  Ffitc2 Ffit
   save sigmafitc2 sigmafit
end

disp('free run..')
dim=length(dkfit);
ufree=zeros(dim, N);
uold=zeros(dim,1);
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end    
   u1=uold+(-dkfit+1i*omegakfit).*uold*dt+Ffit*dt+sqrt(dt)*sigmafit*randn(dim,1);
   uold=u1;
   ufree(:,i)=u1;
end
ufree1=ufree([1,2:2:end],:);
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

for i=1:7
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


% return
figure()
for i=1:7
%     
    subplot(7,2,i*2-1)
    [prob,r1]=ksdensity(real(traj(i,range)));
    plot(r1,prob,'k','Linewidth',1.5);
    [prob,r1]=ksdensity(real(ufree1(i,range)));hold on
    plot(r1,prob,'b','Linewidth',1.5);
    setgca(16)
    
   subplot(7,2,i*2)
    if i==1
    acf1=autocorr(real(traj(i,range)),Naut*10);
    plot(0:dt:Naut*10*dt,acf1,'k','Linewidth',1.5);
    acf1=autocorr(real(ufree1(i,range)),Naut*10);hold on
    plot(0:dt:Naut*10*dt,acf1,'b','Linewidth',1.5);    
    else
     acf1=autocorr(real(traj(i,range)),Naut);
    plot(0:dt:Naut*dt,acf1,'k','Linewidth',1.5); 
    acf1=autocorr(real(ufree1(i,range)),Naut);hold on
    plot(0:dt:Naut*dt,acf1,'b','Linewidth',1.5);    
    end    
    setgca(16)

end
  lgnd=legend('Truth',...
    'Free run fit');
set(lgnd,'FontSize',12);
 filenameStrPdf = sprintf('oufitcase%d.pdf',testcase);   
print(gcf, '-dpdf', '-r600', filenameStrPdf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gapp=100;trajg=traj(:,gapp:gapp:end);Dt=gapp*dt;
figure()
for i=1:3
%     figure()
    subplot(3,3,i*3-2)
    plot(Dt:Dt:Dt*length(trajg),real(trajg(i,:)),'k','linewidth',1);
    setgca(16)
    
    subplot(3,3,i*3-1)
    [prob,r1]=ksdensity(real(traj(i,range)));
    plot(r1,prob,'k','Linewidth',1.5);
    [prob,r1]=ksdensity(real(ufree1(i,range)));hold on
    plot(r1,prob,'b','Linewidth',1.5);
    setgca(16)
    
   subplot(3,3,i*3)
    if i==1
    acf1=autocorr(real(traj(i,range)),Naut*10);
    plot(0:dt:Naut*10*dt,acf1,'k','Linewidth',1.5);
    acf1=autocorr(real(ufree1(i,range)),Naut*10);hold on
    plot(0:dt:Naut*10*dt,acf1,'b','Linewidth',1.5);    
    else
     acf1=autocorr(real(traj(i,range)),Naut);
    plot(0:dt:Naut*dt,acf1,'k','Linewidth',1.5); 
    acf1=autocorr(real(ufree1(i,range)),Naut);hold on
    plot(0:dt:Naut*dt,acf1,'b','Linewidth',1.5);    
    end    
    setgca(16)

end
  lgnd=legend('Truth',...
    'Free run fit');
set(lgnd,'FontSize',12);
 filenameStrPdf = sprintf('trajcase%d.pdf',testcase);   
print(gcf, '-dpdf', '-r600', filenameStrPdf)




return







id=2;ts=1;dt1=0.01;
plotacf1(real(traj(id,ts/dt1:end)),10000,dt1);

% load a;figure();plot(a(2,:),a(1,:),'b')

fitfun=0;
dka=zeros(dim,1);omegaka=zeros(dim,1);mua=zeros(dim,1);sigmak2a=zeros(dim,1);
for i=1:size(traj,1)
    if i==1;
        maxlag1=50000;
    else
        maxlag1=5000;
    end
[dka(i),omegaka(i),mua(i),sigmak2a(i)]=fit_ou_single(funfit,maxlag1,dt,ufree(i,ts/dt:end));  
    
end

load s;acf1=plotacf1(real(s),40000,dt1);
