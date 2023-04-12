function [prob,r1]=plotpdf0(data)
data=data(:);
 [prob,r1]=ksdensity(data);
plot(r1,prob,'k','Linewidth',2); set(gca,'fontsize',16)