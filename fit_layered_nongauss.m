% clear;
close all
% load traj; 

    load dkfitc2;load omegakfitc2;load Ffitc2;load sigmafitc2
dkfit=dkfit([1,2:2:end]);omegakfit=omegakfit([1,2:2:end]);Ffit=Ffit([1,2:2:end]);
sigmafit=sigmafit([1,2:2:end],[1,2:2:end]);

    rng(11)
    ns=99; 
    
  for   ind=1:7
      if ind==1
          traj(1,:)=real(traj(1,:))+1i*real(traj(1,:));
      end
      fprintf('mode is %d\n',ind-1);
lambda=dkfit(ind);%%%%%%%%%%%%% select modes for  

v_truth=real(traj(ind,:));


N=length(v_truth);dt=1/1000;
m = mean(v_truth);%%% 

[fi_v,xx_v] = ksdensity(v_truth);
mi=min(xx_v);ma=max(xx_v);xi=mi: (ma-mi)/ns:ma;
[fi_v,xx_v] = ksdensity(v_truth,xi);
% return
fi_v(fi_v<=1e-5) = 1e-5;
fi_v = fi_v/trapz(xx_v,fi_v);
Delta_v = xx_v(2)-xx_v(1);
Phi = zeros(size(xx_v));
for i = 2:length(Phi)
    Phi(i) = Phi(i-1) + ( (xx_v(i-1) - m) * fi_v(i-1) +...
        (xx_v(i) - m) * fi_v(i) )/2 * (xx_v(i) - xx_v(i-1));
end
% return
sigma_Z_temp = -2 * lambda * Phi./fi_v;
sigma_Z_temp(sigma_Z_temp<0) = 0;
sigma_Z = sqrt(sigma_Z_temp);

        sigma_wx_var_real(ind,:)=sigma_Z;
    xxvx_real(ind,:)=xx_v(1);
    Delta_vx_real(ind)=Delta_v;

%%%%%%
v_truth=imag(traj(ind,:));        

m = mean(v_truth);%%% 

[fi_v_imag,xx_v_imag] = ksdensity(v_truth);
mi=min(xx_v_imag);ma=max(xx_v_imag);xi=mi: (ma-mi)/ns:ma;
[fi_v_imag,xx_v_imag] = ksdensity(v_truth,xi);

fi_v_imag(fi_v_imag<=1e-5) = 1e-5;
fi_v_imag = fi_v_imag/trapz(xx_v_imag,fi_v_imag);

% figure();plot(xx_v_imag,fi_v_imag)

Delta_v_imag = xx_v_imag(2)-xx_v_imag(1);
Phi_imag = zeros(size(xx_v_imag));
for i = 2:length(Phi_imag)
    Phi_imag(i) = Phi_imag(i-1) + ( (xx_v_imag(i-1) - m) * fi_v_imag(i-1) + ...
        (xx_v_imag(i) - m) * fi_v_imag(i) )/2 * (xx_v_imag(i) - xx_v_imag(i-1));
end

sigma_Z_temp = -2 * lambda * Phi_imag./fi_v_imag;

% return
sigma_Z_temp(sigma_Z_temp<0) = 0;
sigma_Z_imag = sqrt(sigma_Z_temp);

    sigma_wx_var_imag(ind,:)=sigma_Z_imag;
    xxvx_imag(ind,:)=xx_v_imag(1);
    Delta_vx_imag(ind)=Delta_v_imag;
  end
  
%   return
  
Z = zeros(1,N*1);
Z1=zeros(1,N*1);
dt1=dt;
for i = 2:N
    sigma_Z_index = round( (real(Z(i-1)) - xx_v(1)) /Delta_v);
    if sigma_Z_index <= 0 || sigma_Z_index >= length(Phi)
        sigma_Z_i = 0;
    else
        sigma_Z_i = sigma_Z(sigma_Z_index);
    end
    
        sigma_Z_index = round( (imag(Z(i-1)) - xx_v_imag(1)) /Delta_v_imag);
    if sigma_Z_index <= 0 || sigma_Z_index >= length(Phi_imag)
        sigma_Z_i_imag = 0;
    else
        sigma_Z_i_imag = sigma_Z_imag(sigma_Z_index);
    end
%     a=sqrt(sigma_Z_i^2+sigma_Z_i_imag^2);
if ind~=1
    Z(i) = Z(i-1) + (-lambda+sqrt(-1)*omegakfit(ind)) * (Z(i-1)) * dt1 + 1*dt1*Ffit(ind)+ ...
        (sigma_Z_i* randn+sqrt(-1)*sigma_Z_i_imag* randn) * sqrt(dt1) ;
else
     Z(i) = Z(i-1) + (-lambda+sqrt(-1)*omegakfit(ind)) * (Z(i-1)) * dt1 + 1*dt1*Ffit(ind)+ ...
        (sigma_Z_i* randn) * sqrt(dt1) ;
   
end
    if ind~=1

    Z1(i) = Z1(i-1) + (-lambda+sqrt(-1)*omegakfit(ind)) * (Z1(i-1)) * dt1 + 1*dt1*Ffit(ind)+ ...
        (sigmafit(ind,ind)* randn+sqrt(-1)*sigmafit(ind,ind)* randn) * sqrt(dt1) ;

    else
    Z1(i) = Z1(i-1) + (-lambda+sqrt(-1)*omegakfit(ind)) * (Z1(i-1)) * dt1 + 1*dt1*Ffit(ind)+ ...
        (sigmafit(ind,ind)* randn) * sqrt(dt1) ;
        
    end
% 
end

% save Delta_vx_real  Delta_vx_real
% save sigma_wx_var_real sigma_wx_var_real
% save xxvx_real xxvx_real
% save Delta_vx_imag Delta_vx_imag
% save sigma_wx_var_imag sigma_wx_var_imag
% save xxvx_imag xxvx_imag
% return
v_truth=real(traj(ind,:));

[prob,r]=ksdensity(v_truth);
[prob1,r1]=ksdensity(real(Z));
[prob2,r2]=ksdensity(real(Z1));

a=[prob;r];save a a
b=[prob1;r1];save b b
c=[prob2;r2];save c c

load a;load b;load c
figure();
subplot(2,1,1)
plot(a(2,:),a(1,:),'k','linewidth',1.5);hold on
plot(b(2,:),b(1,:),'r','linewidth',1.5);hold on
plot(c(2,:),c(1,:),'b','linewidth',1.5);hold on
 title([' K is ',num2str(ind-1) ],'fontsize',14);
setgca(16)

if ind~=1
v_truth=imag(traj(ind,:));
[prob,r]=ksdensity(v_truth);
[prob1,r1]=ksdensity(imag(Z));
[prob2,r2]=ksdensity(imag(Z1));

a1=[prob;r];save a1 a1
b1=[prob1;r1];save b1 b1
c1=[prob2;r2];save c1 c1



load a1;load b1;load c1
subplot(2,1,2)
plot(a1(2,:),a1(1,:),'k','linewidth',1.5);hold on
plot(b1(2,:),b1(1,:),'r','linewidth',1.5);hold on
plot(c1(2,:),c1(1,:),'b','linewidth',1.5);hold on
end
  lgnd=legend('Truth',...
    'non gauss','gauss');
set(lgnd,'FontSize',12);
setgca(16)
%  title([' K is ',num2str(ind-1) ],'fontsize',16);
 filenameStrPdf = sprintf('nongausspdfmode%d.pdf',ind-1);   
print(gcf, '-dpdf', '-r600', filenameStrPdf)






figure()
Naut=2000;if ind==1;Naut=Naut*200;end
subplot(2,1,1)
acf1=autocorr(real(traj(ind,:)),Naut);
 plot(0:dt:Naut*dt,acf1,'k','Linewidth',1.5); hold on
acf1=autocorr(real(Z),Naut);
 plot(0:dt:Naut*dt,acf1,'r','Linewidth',1.5); hold on
 acf1=autocorr(real(Z1),Naut);
 plot(0:dt:Naut*dt,acf1,'b','Linewidth',1.5); setgca(16);
  title([' K is ',num2str(ind-1) ],'fontsize',14);
 if ind~=1
subplot(2,1,2)
acf1=autocorr(imag(traj(ind,:)),Naut);
 plot(0:dt:Naut*dt,acf1,'k','Linewidth',1.5); hold on
acf1=autocorr(imag(Z),Naut);
 plot(0:dt:Naut*dt,acf1,'r','Linewidth',1.5); hold on
 acf1=autocorr(imag(Z1),Naut);
 plot(0:dt:Naut*dt,acf1,'b','Linewidth',1.5); setgca(16);
setgca(16)
 end
   lgnd=legend('Truth',...
    'non gauss','gauss');
set(lgnd,'FontSize',12);


  filenameStrPdf = sprintf('nongaussacfmode%d.pdf',ind-1);   
print(gcf, '-dpdf', '-r600', filenameStrPdf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gapp=100;trajg=traj(ind,gapp:gapp:end);Dt=gapp*dt;

figure();subplot(3,1,1)
plot(Dt:Dt:Dt*length(trajg),real(trajg),'k','linewidth',1);hold on
 title([' K is ',num2str(ind-1) ],'fontsize',14);
subplot(3,1,2);plot(Dt:Dt:Dt*length(trajg),real(Z(gapp:gapp:end)),'r','linewidth',1);hold on
subplot(3,1,3);plot(Dt:Dt:Dt*length(trajg),real(Z1(gapp:gapp:end)),'b','linewidth',1);

%    lgnd=legend('Truth',...
%     'non gauss','gauss');
% set(lgnd,'FontSize',12);
  filenameStrPdf = sprintf('nongausstrajmode%d.pdf',ind-1);   
print(gcf, '-dpdf', '-r600', filenameStrPdf)
