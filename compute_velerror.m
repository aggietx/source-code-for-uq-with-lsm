% function e=compute_velerror(ts,
datacase=3;ii=6;
if datacase==1
u_hatapp=gamma_mean_trace;
% u_hatapp=all_filmean(:,:,ii);
elseif datacase==2
u_hatapp=gamma_mean_trace_smoother;
%  u_hatapp=all_smmean(:,:,ii);
else
u_hatapp=gamma_mean_trace_smoother_sample;
%  u_hatapp=all_bsmean(:,:,ii);
end

ts=5;
rmsvel=zeros(T-ts+1,1);
ccvel=rmsvel;
for time=ts:T
    
       u = (Ga*  u_hat(:,time/dt)); 
    v = (Gb*  u_hat(:,time/dt)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity, record exact'  )
      max(abs(imag(u)))
    end
    u=real(u);v=real(v);
 ue=reshape(u, Dim_Grid,Dim_Grid);
 ve=reshape(v, Dim_Grid,Dim_Grid); 
vex=sqrt(ue.^2+ve.^2);

        u = (Ga*  u_hatapp(:,time/dt)); 
    v = (Gb*  u_hatapp(:,time/dt)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity, record exact'  )
      max(abs(imag(u)))
    end
    u=real(u);v=real(v);
 uflowxapp=reshape(u, Dim_Grid,Dim_Grid);
 vflowyapp=reshape(v, Dim_Grid,Dim_Grid); 
 
      vapp=sqrt(uflowxapp.^2+vflowyapp.^2);
rmsvel(time-ts+1)=rmse(vapp(:),vex(:),0)/std(vex(:));
ccvel(time-ts+1)=corr2(vapp(:),vex(:));
end
fprintf('rms & cc are %2.3f %2.3f\n',mean(rmsvel),mean(ccvel));
figure();subplot(2,1,1);
plot(ts:T,rmsvel,'*-g','linewidth',2);title('RMS of velocity','fontsize',16);setgca(18);xlim([ts,T])
subplot(2,1,2);
plot(ts:T,ccvel,'*-g','linewidth',2);title('CC of velocity','fontsize',16);setgca(18);xlim([ts,T])
xlabel('t','fontsize',20)



