% plotacf(y_truth,gamma_mean_trace(1,:),gamma_mean_trace_smoother(1,:),gamma_mean_trace_smoother_sample(1,:),1:N,Tacf,dt)
close all;
Linewid=2;
%%% L63
% figure()
% subplot(3,1,1)
% p1 = plot(dt:dt:T,x_truth,'red','Linewidth',Linewid);
% xlim([0,50])
% title('Ture signal of x','FontSize',18)
% subplot(3,1,2)
% p2 = plot(dt:dt:T,y_truth,'r');
% hold on;plot(dt:dt:T,gamma_mean_trace(1,:),'magenta','Linewidth',Linewid);
% plot(dt:dt:T,gamma_mean_trace_smoother(1,:),'blue','Linewidth',Linewid);
% plot(dt:dt:T,gamma_mean_trace_smoother_sample(1,:),'green','Linewidth',Linewid);
% xlim([0,50])
% title('Trajectory comparison of y','FontSize',18)
% subplot(3,1,3)
% p3 = plot(dt:dt:T,z_truth,'r');
% hold on;plot(dt:dt:T,gamma_mean_trace(2,:),'magenta','Linewidth',Linewid);
% plot(dt:dt:T,gamma_mean_trace_smoother(2,:),'blue','Linewidth',Linewid);
% plot(dt:dt:T,gamma_mean_trace_smoother_sample(2,:),'green','Linewidth',Linewid);
% xlim([0,50])
% title('Trajectory comparison of z','FontSize',18)
% lgd = legend({'Truth','Filtered mean','Smoothed mean','Backward sampling'},...
%   'Orientation','horizontal',  'FontSize',18,'Location', 'southeast');
% xlabel('t', 'FontSize',18)
% return
%%%% L84
%%%% L84
%%%% L84
ind=21;
figure()
ranget=50/dt:100/dt;

subplot(3,2,1)
p1 = plot(50:dt:100,u_truth(ind,ranget),'red','Linewidth',Linewid);
box on
xlim([50,100])
title('Ture signal of u_{21}','FontSize',18)
% xlabel('t', 'FontSize',18)

subplot(3,2,2)
p2 = plot(50:dt:100,v_truth(ind,ranget),'r');
hold on;plot(50:dt:100,gamma_mean_trace(ind,ranget),'magenta','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother(ind,ranget),'blue','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother_sample(ind,ranget),'green','Linewidth',Linewid);
box on
xlim([50,100])
ylim([-11,-5])
title('Trajectory comparison of v_{21,1} ','FontSize',18)
% xlabel('t', 'FontSize',18)
ind=21+40;
subplot(3,2,3)
p3 = plot(50:dt:100,v_truth(ind,ranget),'r','Linewidth',Linewid);
hold on;plot(50:dt:100,gamma_mean_trace(ind,ranget),'magenta','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother(ind,ranget),'blue','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother_sample(ind,ranget),'green','Linewidth',Linewid);
box on
xlim([50,100])
title('Trajectory comparison of v_{21,2} ','FontSize',18)
% xlabel('t', 'FontSize',18)
ind=21+80;
subplot(3,2,4)
p4 = plot(50:dt:100,v_truth(ind,ranget),'r');
hold on;plot(50:dt:100,gamma_mean_trace(ind,ranget),'magenta','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother(ind,ranget),'blue','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother_sample(ind,ranget),'green','Linewidth',Linewid);
box on
xlim([50,100])
title('Trajectory comparison of v_{21,3} ','FontSize',18)
% xlabel('t', 'FontSize',18)
ind=21+120;
subplot(3,2,5)
p5= plot(50:dt:100,v_truth(ind,ranget),'r');
hold on;plot(50:dt:100,gamma_mean_trace(ind,ranget),'magenta','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother(ind,ranget),'blue','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother_sample(ind,ranget),'green','Linewidth',Linewid);
box on
xlim([50,100])
title('Trajectory comparison of v_{21,4} ','FontSize',18)
xlabel('t', 'FontSize',18)
ind=21+160;
subplot(3,2,6)
p6 = plot(50:dt:100,v_truth(ind,ranget),'r');
hold on;plot(50:dt:100,gamma_mean_trace(ind,ranget),'magenta','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother(ind,ranget),'blue','Linewidth',Linewid);
plot(50:dt:100,gamma_mean_trace_smoother_sample(ind,ranget),'green','Linewidth',Linewid);
box on
xlim([50,100])
title('Trajectory comparison of v_{21,5} ','FontSize',18)
xlabel('t', 'FontSize',18)
lgd = legend({'Truth','Filtered mean','Smoothed mean','Backward sampling'},...
  'Orientation','horizontal',  'FontSize',12,'Location', 'southeast');
