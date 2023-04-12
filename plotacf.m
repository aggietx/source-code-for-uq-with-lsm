function plotacf(v_truth,filter_mean,smoothed_mean,back_sample,range,t,dt)
figure();[r1]=autocorr(v_truth(range),floor(t/dt));
[r2]=autocorr(back_sample(range),floor(t/dt));
[r3]=autocorr(smoothed_mean(range),floor(t/dt));
[r4]=autocorr(filter_mean(range),floor(t/dt));
plot(0:dt:dt*(length(r1)-1),r1,'r','Linewidth',2);hold on;
plot(0:dt:dt*(length(r4)-1),r4,'magenta','Linewidth',2)
plot(0:dt:dt*(length(r3)-1),r3,'blue','Linewidth',2)
plot(0:dt:dt*(length(r2)-1),r2,'green','Linewidth',2)
legend({'Truth','Filtered mean','Smoothed mean','Backward sampling'},'FontSize',12)
xlabel('t')

