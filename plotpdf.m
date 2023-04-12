function [prob,r1]=plotpdf(data)
data=data(:);
 [prob,r1]=ksdensity(data);
figure();plot(r1,prob,'r','Linewidth',2); set(gca,'fontsize',18)
title('PDF','fontsize',18)