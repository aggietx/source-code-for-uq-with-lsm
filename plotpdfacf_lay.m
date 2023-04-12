clear;
% close all
load trajc2;load ufreec2
% load trajc1;load ufreec1
dt=0.01;ts=10;Naut=1000;
range=ts/dt:5000/dt;
figure()
for i=1:7
    subplot(7,2,2*i-1)
    [prob,r1]=ksdensity(real(traj(i,range)));
    plot(r1,prob,'k','Linewidth',2);
    [prob,r1]=ksdensity(real(ufree(i,range)));hold on
    plot(r1,prob,'b','Linewidth',2);
    setgca(16)
    
    subplot(7,2,2*i)
    if i==1
    acf1=autocorr(real(traj(i,range)),Naut*10);
    plot(0:dt:Naut*10*dt,acf1,'k','Linewidth',2);
    acf1=autocorr(real(ufree(i,range)),Naut*10);hold on
    plot(0:dt:Naut*10*dt,acf1,'b','Linewidth',2);    
    else
     acf1=autocorr(real(traj(i,range)),Naut);
    plot(0:dt:Naut*dt,acf1,'k','Linewidth',2); 
    acf1=autocorr(real(ufree(i,range)),Naut);hold on
    plot(0:dt:Naut*dt,acf1,'b','Linewidth',2);    
    end    
    setgca(16)
end
  lgnd=legend('Truth',...
    'Free run fit');
set(lgnd,'FontSize',12);