subplot(6,1,2)
p2 = plot(dt:dt:T,v_truth(ind,:),'r');
hold on;plot(dt:dt:T,gamma_mean_trace(ind,:),'magenta');
plot(dt:dt:T,gamma_mean_trace_smoother(ind,:),'blue');
plot(dt:dt:T,gamma_mean_trace_smoother_sample(ind,:),'green');
box on
xlim([50,100])
title('Trajectory comparison of v_{21,1} ','FontSize',6)

ind=21+40;
subplot(6,1,3)
p3 = plot(dt:dt:T,v_truth(ind,:),'r');
hold on;plot(dt:dt:T,gamma_mean_trace(ind,:),'magenta');
plot(dt:dt:T,gamma_mean_trace_smoother(ind,:),'blue');
plot(dt:dt:T,gamma_mean_trace_smoother_sample(ind,:),'green');
box on
xlim([50,100])
title('Trajectory comparison of v_{21,2} ','FontSize',6)

ind=21+80;
subplot(6,1,4)
p4 = plot(dt:dt:T,v_truth(ind,:),'r');
hold on;plot(dt:dt:T,gamma_mean_trace(ind,:),'magenta');
plot(dt:dt:T,gamma_mean_trace_smoother(ind,:),'blue');
plot(dt:dt:T,gamma_mean_trace_smoother_sample(ind,:),'green');
box on
xlim([50,100])
title('Trajectory comparison of v_{21,3} ','FontSize',6)

ind=21+120;
subplot(6,1,5)
p5= plot(dt:dt:T,v_truth(ind,:),'r');
hold on;plot(dt:dt:T,gamma_mean_trace(ind,:),'magenta');
plot(dt:dt:T,gamma_mean_trace_smoother(ind,:),'blue');
plot(dt:dt:T,gamma_mean_trace_smoother_sample(ind,:),'green');
box on
xlim([50,100])
title('Trajectory comparison of v_{21,4} ','FontSize',6)

ind=21+160;
subplot(6,1,6)
p6 = plot(dt:dt:T,v_truth(ind,:),'r');
hold on;plot(dt:dt:T,gamma_mean_trace(ind,:),'magenta');
plot(dt:dt:T,gamma_mean_trace_smoother(ind,:),'blue');
plot(dt:dt:T,gamma_mean_trace_smoother_sample(ind,:),'green');
box on
xlim([50,100])
title('Trajectory comparison of v_{21,5} ','FontSize',6)
xlabel('t')
lgd = legend({'Truth','Filtered mean','Smoothed mean','Backward sampling'},...
  'Orientation','horizontal',  'FontSize',6,'Location', 'southeast');
