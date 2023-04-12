ts=5;gap=50;
ix=3;iy=3;ind=getind(ix,iy,kk);
figure();plot(ts:dt*gap:T,real(u_hat(ind,ts/dt:gap:end)),'r','linewidth',1);hold on;
plot(ts:dt*gap:T,real(gamma_mean_trace(ind,ts/dt:gap:end)),'b','linewidth',1);
plot(ts:dt*gap:T,real(gamma_mean_trace_smoother(ind,ts/dt:gap:end)),'g','linewidth',1);
plot(ts:dt*gap:T,real(gamma_mean_trace_smoother_sample(ind,ts/dt:gap:end)),'k','linewidth',1);
lgnd=legend('exact','filter','smoother','sampled');
set(lgnd,'FontSize',16);
setgca(16);title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
xlabel('t','fontsize',20)
xlim([ts,T])