L1=2*K_max+1;
F=fft(eye(2*K_max+1));
fft2mat = kron(F,F);
Dim_Grid=2*K_max+1;
xy=kk0'/(2*K_max+1)*2*pi;
data=u_hat(:,100);
rketa=[1./sqrt(kk0(1,:).^2 + kk0(2,:).^2+1),...
    1./sqrt(kk0(1,:).^2 + kk0(2,:).^2)/sqrt(2)./sqrt(kk0(1,:).^2 + kk0(2,:).^2 + 1).*(kk0(1,:).^2 + kk0(2,:).^2),...
    1./sqrt(kk0(1,:).^2 + kk0(2,:).^2)/sqrt(2)./sqrt(kk0(1,:).^2 + kk0(2,:).^2 + 1).*(kk0(1,:).^2 + kk0(2,:).^2)];
rketa(1+dimuhat0)=0;rketa(1+dimuhat0*2)=0;
    invG1 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rk(1,:))); % Fourier bases for u
    invG2 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rk(2,:))); % Fourier bases for v
    invG3 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rketa)); % Fourier bases for v
    u = (invG1* data); 
    v = (invG2* data); 
    eta=(invG3* data);

  G=exp(-1i  * xy*kk0)/(2*K_max+1)^2;%%%%% shallow water matrix from physical domain to F domain
  
    data1=zeros(size(data));
  ub=G*u;
  vb=G*v;
  etab=G*eta;
    for ii=2:length(kk0)
       ind=[ii,ii+dimuhat0,ii+2*dimuhat0];
       ux=ub(ii);vx=vb(ii);etax=etab(ii);
       rmatrix=[rk(:,ind);rketa(:,ind)];
       data1(ind)=rmatrix\[ux;vx;etax];
    end
%     data1([1,1+dimuhat0,1+2*dimuhat0])=0;
    norm(data-data1)
return
%% backward sampling
disp('backward sampling...')
time=1;
gamma_mean0_smoother_sample = gamma_mean0;
gamma_mean_trace_smoother_sample = zeros(Dim_Y,N);
gamma_mean_trace_smoother_sample(redind,end)=gamma_mean0;
for i=N:-1:2
   if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',i);
    end    
       R=gamma_cov_trace(:,:,i);R=regularize(R,regularizeval);
       invR=R\speye(Dim_Yr);

    gamma_mean_smoother_sample=gamma_mean0_smoother_sample+dt*(-a0-a1*gamma_mean0_smoother_sample+...
        (b1 * b1')*invR*(gamma_mean_trace_red(:,i)-gamma_mean0_smoother_sample))+sqrt(dt)*b1*randn(Dim_Yr,1);
gamma_mean_trace_smoother_sample(redind,i-1) = gamma_mean_smoother_sample;

    gamma_mean0_smoother_sample = gamma_mean_smoother_sample;
    if mod(i,1/dt)==0
    u = (Ga*  gamma_mean_trace_smoother_sample(:,i-1)); 
    v = (Gb*  gamma_mean_trace_smoother_sample(:,i-1)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity,record backward sampling'  )
      max(abs(imag(u)))
    end
    uexactbackward(:,time)=real(u);
    vexactbackward(:,time)=real(v);
     
        time=time+1;  
    end    
end
figure();plot(dt:dt:T,real(u_hat(ind,:)),'r','linewidth',1.5);hold on;
plot(dt:dt:T,real(gamma_mean_trace_smoother_sample(ind,:)),'b','linewidth',1.5);
setgca(16);title(['Backward sampling, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
xlabel('t','fontsize',16)

rnorm(real(gamma_mean_trace(ind,:)),real(u_hat(ind,:)));
rnorm(real(gamma_mean_trace_smoother(ind,:)),real(u_hat(ind,:)));
rnorm(real(gamma_mean_trace_smoother_sample(ind,:)),real(u_hat(ind,:)));


corr2(real(gamma_mean_trace(ind,:)),real(u_hat(ind,:)));
corr2(real(gamma_mean_trace_smoother(ind,:)),real(u_hat(ind,:)));
corr2(real(gamma_mean_trace_smoother_sample(ind,:)),real(u_hat(ind,:)));