function acf1=plotacf1(data,Naut,dt)

 data=data(:);

acf1 = autocorr(data,Naut);


figure()
plot(0:dt:Naut*dt,acf1,'r','linewidth',2);
set(gca,'fontsize',18)
title('ACF','fontsize',18)
xlabel('t','fontsize',18)