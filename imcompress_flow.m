%%%% filter and smoother for impressible flow
clear;close all;rng(1);
dt=1/500;T=20;N=T/dt;regularizeval=10^(-4);
K_max=3;L=24;
E_0=1;k_0=2;alpha=3;nosmoother=1;
dimuhat=(2*K_max+1)^2;Kmax=K_max;
[ky,kx]=meshgrid([0:K_max,-K_max:-1],[0:K_max,-K_max:-1]);
kk=[kx(:),ky(:)]';
rk=zeros(size(kk));
for i=2:dimuhat
        ki=kk(:,i);
        rk(:,i)=1i*[-ki(2);ki(1)]/norm(ki);
end
ix=Kmax;iy=Kmax;ind=getind(ix,iy,kk);%% sample plot

L1=2*pi;%%% size of the domain
xexact=zeros(L,N);
yexact=xexact;
xexact(:,1)=L1*rand(L,1);
yexact(:,1)=L1*rand(L,1);
sig_ex=0.25;
u_hat=zeros(dimuhat,N);
% % dkexact=-0.5*ones(dimuhat,1);
% % omegaexact0=0.5;omegaexact=formomega(dimuhat,omegaexact0,kk);
% % Fuexact=0.2*ones(dimuhat,1);
% % sigmauexact=0.5*ones(dimuhat,1);
% % sigmauhat=form_sigmatrix1(dimuhat,sigmauexact,kk);

[dkexact,omegaexact,Fuexact,sigmauhat,ene,sigmak2exact]=form_coeff(dimuhat,kk,E_0,k_0,alpha,Kmax);
dkexact=-dkexact;


%% generate true signal
disp('generating true signal...')
for i=2:N
    u_hat(:,i)=u_hat(:,i-1)+(dkexact+1i*omegaexact).*u_hat(:,i-1)*dt+...
        Fuexact*dt+sqrt(dt)*sigmauhat*randn(dimuhat,1);
    if max(abs(real(u_hat(:,i))))>10^8
        disp('error, data blow up')
    end
end
Dim_Grid = 40;
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
uexact=zeros(Dim_Grid^2,T);vexact=zeros(Dim_Grid^2,T);
uexactfilter=zeros(Dim_Grid^2,T);vexactfilter=zeros(Dim_Grid^2,T);
uexactsmoother=zeros(Dim_Grid^2,T);vexactsmoother=zeros(Dim_Grid^2,T);
uexactbackward=zeros(Dim_Grid^2,T);vexactbackward=zeros(Dim_Grid^2,T);
    Ga = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(1,:))); % Fourier bases for u
    Gb = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(2,:))); % Fourier bases for v
% return
time=1;
for i=2:N
    x_loc = [xexact(:,i-1),yexact(:,i-1)];
    G1 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
    G2 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
    u = (G1*  u_hat(:,i-1)); 
    v = (G2*  u_hat(:,i-1)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity'  )
      max(abs(imag(u)))
    end
    u=real(u);
    v=real(v);
    xexact(:,i) = xexact(:,i-1) + u * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in x
    yexact(:,i) = yexact(:,i-1) + v * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in y
%     xexact=real(xexact);yexact=real(yexact);
if imag(max(abs(xexact(:,i))))>0 || imag(max(abs(yexact(:,i))))>0
    disp('error, complex')
end
    xexact(:,i) = mod(xexact(:,i),L1);
    yexact(:,i) = mod(yexact(:,i),L1);
    if mod(i,1/dt)==0
    u = (Ga*  u_hat(:,i-1)); 
    v = (Gb*  u_hat(:,i-1)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity, record exact'  )
      max(abs(imag(u)))
    end
    uexact(:,time)=real(u);
    vexact(:,time)=real(v);
     
        time=time+1;  
    end
end
% return
%% filter
Dim_Y=dimuhat;Dim_X=2*L;
gamma_mean0 = zeros(Dim_Y,1);
gamma_cov0 = eye(Dim_Y)*0.01; 
gamma_mean_trace=zeros(Dim_Y,N);
gamma_mean_trace(:,1) = gamma_mean0;
if nosmoother~=1
gamma_cov_trace=zeros(Dim_Y,Dim_Y,N);
gamma_cov_trace(:,:,1) = gamma_cov0;
end
A0=zeros(Dim_X,1);
a0=Fuexact;
a1=diag(dkexact+1i*omegaexact);
b1=sigmauhat;
x=xexact;y=yexact;
invBoB = 1 /sig_ex/ sig_ex * eye(2*L); % inverse of the square of the observational noise
time=1;
disp('data assimilation (filter)......')
for i = 2:N
    if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
    end
    % observational operator 
    x_loc = [xexact(:,i-1),yexact(:,i-1)];
   
    G1 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
    G2 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
    % computing the difference between the locations in the Lagrangian
    % tracers; need to consider the cases near the boundaries 
    diff_x1 = x(:,i) - x(:,i-1); diff_x2 = x(:,i) - x(:,i-1) + L1; diff_x3 = x(:,i) - x(:,i-1) - L1;  
    diff_y1 = y(:,i) - y(:,i-1); diff_y2 = y(:,i) - y(:,i-1) + L1; diff_y3 = y(:,i) - y(:,i-1) - L1;  
    diff_xtemp = min(abs(diff_x1), abs(diff_x2)); diff_x_index = min(abs(diff_x3), diff_xtemp);
    diff_ytemp = min(abs(diff_y1), abs(diff_y2)); diff_y_index = min(abs(diff_y3), diff_ytemp);
    diff_x1_index = (diff_x_index == abs(diff_x1)); diff_x2_index = (diff_x_index == abs(diff_x2)); diff_x3_index = (diff_x_index == abs(diff_x3)); 
    diff_y1_index = (diff_y_index == abs(diff_y1)); diff_y2_index = (diff_y_index == abs(diff_y2)); diff_y3_index = (diff_y_index == abs(diff_y3)); 
    diff_x = diff_x1 .* diff_x1_index + diff_x2 .* diff_x2_index + diff_x3 .* diff_x3_index;
    diff_y = diff_y1 .* diff_y1_index + diff_y2 .* diff_y2_index + diff_y3 .* diff_y3_index;
    diff_xy = [diff_x; diff_y];
    
    A1=[G1;G2];
    


    % run the data assimilation for posterior mean and posterior covariance
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (diff_xy - A0*dt-A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
%     gamma_mean(1:3*L) = real(gamma_mean(1:3*L));
% gamma_mean(1:12) 
% pause


    % save the posterior statistics
    gamma_mean_trace(:,i) = gamma_mean;
    if nosmoother~=1
    gamma_cov_trace(:,:,i) = gamma_cov;
    end
    % update
    gamma_mean0 = gamma_mean;
    gamma_cov0 = gamma_cov;
%     gamma_cov0 = diag(diag(gamma_cov));
    
    if mod(i,1/dt)==0
    u = (Ga*  gamma_mean); 
    v = (Gb*  gamma_mean); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity,  record filter'  )
      max(abs(imag(u)))
    end
    uexactfilter(:,time)=real(u);
    vexactfilter(:,time)=real(v);
     
        time=time+1;  
    end
end

figure();plot(dt:dt:T,real(u_hat(ind,:)),'r','linewidth',1.5);hold on;
plot(dt:dt:T,real(gamma_mean_trace(ind,:)),'b','linewidth',1.5);
setgca(16);title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
xlabel('t','fontsize',16)
rnorm(real(gamma_mean_trace(ind,:)),real(u_hat(ind,:)));

% imagesc(reshape(uexact(:,5).^2+vexact(:,5).^2,40,40))
if nosmoother
return
end
%% smoother
disp('smoother....')
gamma_mean0_smoother = gamma_mean0;
gamma_cov0_smoother = gamma_cov0;
gamma_mean_trace_smoother = zeros(Dim_Y,N);
gamma_mean_trace_smoother(:,end)=gamma_mean0;

time=1;
for i=N:-1:2
   if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',i);
    end
    % observational operator 
   
    R=gamma_cov_trace(:,:,i);R=regularize(R,regularizeval);
    invR=R\speye(Dim_Y);
 
    gamma_mean_smoother=gamma_mean0_smoother+dt*(-a0-a1*gamma_mean0_smoother+(b1 * b1')*invR*(gamma_mean_trace(:,i)-gamma_mean0_smoother));
    gamma_cov_smoother=gamma_cov0_smoother-dt*( (a1+(b1*b1')*invR)*gamma_cov0_smoother+...
        gamma_cov0_smoother*(a1'+b1*b1'*invR)-b1*b1' );
        
    gamma_mean_trace_smoother(:,i-1) = gamma_mean_smoother;

    gamma_mean0_smoother = gamma_mean_smoother;
    gamma_cov0_smoother = gamma_cov_smoother;

    if mod(i,1/dt)==0
    u = (Ga*  gamma_mean_smoother); 
    v = (Gb*  gamma_mean_smoother); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity, record smoother'  )
      max(abs(imag(u)))
    end
    uexactbackward(:,time)=real(u);
    vexactbackward(:,time)=real(v);
     
        time=time+1;  
    end       
end
figure();plot(dt:dt:T,real(u_hat(ind,:)),'r','linewidth',1.5);hold on;
plot(dt:dt:T,real(gamma_mean_trace_smoother(ind,:)),'b','linewidth',1.5);
setgca(16);title(['Smoother, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
xlabel('t','fontsize',16)
rnorm(real(gamma_mean_trace_smoother(ind,:)),real(u_hat(ind,:)));

%% backward sampling
disp('backward sampling...')
time=1;
gamma_mean0_smoother_sample = gamma_mean0;
gamma_mean_trace_smoother_sample = zeros(Dim_Y,N);
gamma_mean_trace_smoother_sample(:,end)=gamma_mean0;
for i=N:-1:2
   if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',i);
    end    
       R=gamma_cov_trace(:,:,i);R=regularize(R,regularizeval);
       invR=R\speye(Dim_Y);

    gamma_mean_smoother_sample=gamma_mean0_smoother_sample+dt*(-a0-a1*gamma_mean0_smoother_sample+...
        (b1 * b1')*invR*(gamma_mean_trace(:,i)-gamma_mean0_smoother_sample))+sqrt(dt)*b1*randn(Dim_Y,1);
gamma_mean_trace_smoother_sample(:,i-1) = gamma_mean_smoother_sample;

    gamma_mean0_smoother_sample = gamma_mean_smoother_sample;
    if mod(i,1/dt)==0
    u = (Ga*  gamma_mean_smoother_sample); 
    v = (Gb*  gamma_mean_smoother_sample); 
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
