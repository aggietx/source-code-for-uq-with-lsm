clear;
rng(10);
close all;
L=30;K_max=3;L1=2*pi;
dt=0.01;T=100;N=T/dt;
xexact = zeros(L,N);
yexact = zeros(L,N);
% rng(20)
xexact(:,1)=L1*rand(L,1);
% rng(10)
yexact(:,1)=L1*rand(L,1);
sig_ex=0.2;%%%% noise of dx=v*dt
sigma_B = 0.15;0.4/2; % noise of the GB modes
% sigma_g = 0.1;0.4/4; % noise of the gravity modes
% only have GB wave
% 2D shallow water equation
% K_max =1; % the range of Fourier modes is [-K_max, K_max]^2
epsilon = 0.1; % Rossby number
k = zeros(2, (2 * K_max + 1) * (2 * K_max + 1)); % Total number of Fourier wavenumbers
% T = 4; % total time
% if min(thickness)<0.2
% dt = 0.00005; % time step
% else
% dt = 0.0005; % time step
% end
% dt=0.001;
% if K_max == 1;
% dt=0.001;
% end
% arranging Fourier wavenumbers
% arranging in such a way that the complex conjugates modes are next to
% each other, namely (-k1,-k2) will be next to (k1,k2). This will
% facilitate data assimilation and allows it to be a nearly block diagonal
% matrix where each block is 2 by 2.


m = 1;
for i = - K_max : K_max
    if i < 0
        for j = - K_max: i
            k(1, m) = i;
            k(2, m) = j;
            m = m + 2;
        end
    else
        for j = - K_max : i - 1
            k(1, m) = i;
            k(2, m) = j;
            m = m + 2;
        end
    end
end
k(:, 2: 2: end - 1) = - k(:, 1 : 2 : end - 2);
fprintf('kmax is %d\n',K_max);
 kk=k;
[ky,kx]=meshgrid([0:K_max,-K_max:-1],[0:K_max,-K_max:-1]);
kkfft=[kx(:),ky(:)]';
indkk1=zeros( (K_max+1)^2,1);indkk2=zeros((K_max+1)^2,1);
for i=1:(2*K_max+1)^2
indkk1(i)=find(ismember(kkfft',kk(:,i)','rows'));%%%% fft2
indkk2(i)=find(ismember(kk',kkfft(:,i)','rows'));%%%% ifft2 
end
% rk1: GB; 
rk1 = [1./sqrt(k(1,:).^2 + k(2,:).^2+1) .* (-1i * k(2,:));
    1./sqrt(k(1,:).^2 + k(2,:).^2+1) .* (1i * k(1,:))];
rk=rk1;

% stochastic systems for each Fourier mode
% written in the vector form
% dU = (a1 U + a0 ) dt + b1 dW
rng(500)
N = round(T/dt);
Dim_U = length(kk(1,:)); % dimension of the system
Dim_Ug=0;
Dim_UB = Dim_U ;
u_hat = zeros(Dim_U,N); % define all the Fourier modes
% return
d_B = 0.5; % damping of the GB modes
d_g = 0.5; % damping of the gravity modes
% sigma_B = 0.1;0.4/2; % noise of the GB modes
% sigma_g = 0.2;0.4/4; % noise of the gravity modes
f_amp = 0.1;0.005; % large-scale forcing amplitude
f_phase =2*pi/14; 0.5; % large-scale forcing period
f_x_b = 0.5; % large-scale forcing background in x direction
f_y_b = 0;% large-scale forcing background in y direction
% b1: noise coefficient; a1: damping and phase 
b1 = zeros(Dim_U, Dim_U);
b2 = zeros(Dim_U, Dim_U);
a1 = - diag([d_g * ones(1,Dim_Ug * 2), d_g * ones(1,Dim_UB)]) ;% 1i * diag([omegak, zeros(1,Dim_UB+2)]);
L_u = a1;
for i = 1: 2: 2 * Dim_Ug + Dim_UB-1
i1=find(ismember(kkfft',kk(:,i)','rows'));%%%% fft2
i2=find(ismember(kkfft',kk(:,i+1)','rows'));%%%% fft2
        b2(i1,i1) = 1 / sqrt(2) * sigma_B;
        b2(i2, i2) = -1i / sqrt(2) * sigma_B;
        b2(i1, i2) = 1i / sqrt(2) * sigma_B;
        b2(i2, i1) = 1 / sqrt(2) * sigma_B;
        
        b1(i,i) = 1 / sqrt(2) * sigma_B;
        b1(i+1, i+1) = -1i / sqrt(2) * sigma_B;
        b1(i, i+1) = 1i / sqrt(2) * sigma_B;
        b1(i+1, i) = 1 / sqrt(2) * sigma_B;
   
end
Sigma_u = b1;%%% non standard ordering
% Sigma_u= b2;rk=rk(:,indkk2);%%%% standard fft ordering
% numerical integral
rd = zeros(Dim_U,N);
% % rd(1:2:end-2,:) = randn(Dim_Ug + (Dim_UB-1)/2, N) + 1i * randn(Dim_Ug + (Dim_UB-1)/2, N);
% % rd(2:2:end-1,:) = conj(rd(1:2:end-2,:)); % noise needs to be complex conjudate within each 2x2 block
rd(1:end, :) = randn(Dim_UB, N); 
% return
% % % [index1]=find(abs(kk(1,:))<2);[index2]=find(abs(kk(2,:))<2);
% % %  small_index1=intersect(index1,index2);
% % %  small_index1=small_index1(small_index1>2 * Dim_Ug);
% % %  small_index=small_index1(:,1:end-2);
F_u= zeros(Dim_U, 1);
for i = 2:N
%         if mod(i,floor(N/5)) == 0
% %         disp(i*dt)
%  fprintf('rest step is %d\n',N-i);
%         end
    t = i*dt;
   F_u(1:2:end-2) = f_amp * exp(1i * f_phase * t) * ones(floor(Dim_UB/2), 1); 0.4+0.4*1i;%?
   F_u(2:2:end-1) = f_amp * exp(- 1i * f_phase * t) * ones(floor(Dim_UB/2), 1); 0.4-0.4*1i;%?    
F_u(1:end-1)=.1; 
%     a0(end) = 0;f_amp * sin(f_phase * t) + f_y_b; 0;  %%%% should also change in data assimilation code 
    u_hat(:,i) = u_hat(:,i-1) + (L_u * u_hat(:,i-1) + F_u) * dt + ...
        Sigma_u * sqrt(dt) * rd(:, i);
end


for i=2:N
    x_loc = [xexact(:,i-1),yexact(:,i-1)];
    G1 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
    G2 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
    u = real(G1*  u_hat(:,i-1)); 
    v = real(G2*  u_hat(:,i-1)); 
    xexact(:,i) = xexact(:,i-1) + u * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in x
    yexact(:,i) = yexact(:,i-1) + v * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in y
%     xexact=real(xexact);yexact=real(yexact);
if imag(max(abs(xexact(:,i))))>0 || imag(max(abs(yexact(:,i))))>0
    disp('error, complex')
end
    xexact(:,i) = mod(xexact(:,i),L1);
    yexact(:,i) = mod(yexact(:,i),L1);

end

% return
%% filter
Dim_Y=size(u_hat,1);Dim_X=2*L;
gamma_mean0 = zeros(Dim_Y,1);
gamma_cov0 = eye(Dim_Y)*0.01; 
gamma_mean_trace=zeros(Dim_Y,N);gamma_cov_trace=zeros(Dim_Y,Dim_Y,N);
gamma_mean_trace(:,1) = gamma_mean0;
gamma_cov_trace(:,:,1) = gamma_cov0;
A0=zeros(Dim_X,1);
% a0=F_u;
a1=L_u;
x=xexact;y=yexact;
invBoB = 1 /sig_ex/ sig_ex * eye(2*L); % inverse of the square of the observational noise

disp('data assimilation......')
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
    
        t = i*dt;
   F_u(1:2:end-2) = f_amp * exp(1i * f_phase * t) * ones(floor(Dim_UB/2), 1); 0.4+0.4*1i;%?
   F_u(2:2:end-1) = f_amp * exp(- 1i * f_phase * t) * ones(floor(Dim_UB/2), 1); 0.4-0.4*1i;%?    
   F_u(1:end-1)=.1;  
   a0=F_u;
    % run the data assimilation for posterior mean and posterior covariance
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (diff_xy - A0*dt-A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
%     gamma_mean(1:3*L) = real(gamma_mean(1:3*L));
% gamma_mean(1:12) 
% pause


    % save the posterior statistics
    gamma_mean_trace(:,i) = gamma_mean;
    gamma_cov_trace(:,:,i) = gamma_cov;

    % update
    gamma_mean0 = gamma_mean;
%     gamma_cov0 = gamma_cov;
    gamma_cov0 = diag(diag(gamma_cov));
end
plot2line(real(u_hat(12,:)),real(gamma_mean_trace(12,:)));

%% smoother
disp('smoother....')
gamma_mean0_smoother = gamma_mean0;
gamma_cov0_smoother = gamma_cov0;
gamma_mean_trace_smoother = zeros(Dim_Y,N);
gamma_mean_trace_smoother(:,end)=gamma_mean0;

% return
for i=N:-1:2
   if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',i);
    end
    % observational operator 
   
    R=gamma_cov_trace(:,:,i);invR=R\speye(Dim_Y);


            t = i*dt;
   F_u(1:2:end-2) = f_amp * exp(1i * f_phase * t) * ones(floor(Dim_UB/2), 1); 0.4+0.4*1i;%?
   F_u(2:2:end-1) = f_amp * exp(- 1i * f_phase * t) * ones(floor(Dim_UB/2), 1); 0.4-0.4*1i;%?    
  F_u(1:end-1)=.1; 
   a0=F_u;
    
    gamma_mean_smoother=gamma_mean0_smoother+dt*(-a0-a1*gamma_mean0_smoother+(b1 * b1')*invR*(gamma_mean_trace(:,i)-gamma_mean0_smoother));
    gamma_cov_smoother=gamma_cov0_smoother-dt*( (a1+(b1*b1')*invR)*gamma_cov0_smoother+...
        gamma_cov0_smoother*(a1'+b1*b1'*invR)-b1*b1' );
        
    gamma_mean_trace_smoother(:,i-1) = gamma_mean_smoother;

    gamma_mean0_smoother = gamma_mean_smoother;
    gamma_cov0_smoother = gamma_cov_smoother;
end
plot2line(real(u_hat(12,:)),real(gamma_mean_trace_smoother(12,:)));

%% backward sampling
disp('backward sampling...')
gamma_mean0_smoother_sample = gamma_mean0;
gamma_mean_trace_smoother_sample = zeros(Dim_Y,N);
gamma_mean_trace_smoother_sample(:,end)=gamma_mean0;
for i=N:-1:2
   if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',i);
    end    
       R=gamma_cov_trace(:,:,i);invR=R\speye(Dim_Y);
            t = i*dt;
   F_u(1:2:end-2) = f_amp * exp(1i * f_phase * t) * ones(floor(Dim_UB/2), 1); 0.4+0.4*1i;%?
   F_u(2:2:end-1) = f_amp * exp(- 1i * f_phase * t) * ones(floor(Dim_UB/2), 1); 0.4-0.4*1i;%?    
   F_u(1:end-1)=.1; 
   a0=F_u;
    gamma_mean_smoother_sample=gamma_mean0_smoother_sample+dt*(-a0-a1*gamma_mean0_smoother_sample+...
        (b1 * b1')*invR*(gamma_mean_trace(:,i)-gamma_mean0_smoother_sample))+sqrt(dt)*b1*randn(Dim_Y,1);
gamma_mean_trace_smoother_sample(:,i-1) = gamma_mean_smoother_sample;

    gamma_mean0_smoother_sample = gamma_mean_smoother_sample;
end
plot2line(real(u_hat(12,:)),real(gamma_mean_trace_smoother_sample(12,:)));
% plot2line(real(u_hat(1,:)),real(gamma_mean_trace_smoother_sample(1,:)));
s=12;
rnorm(real(gamma_mean_trace(s,:)),real(u_hat(s,:)));
rnorm(real(gamma_mean_trace_smoother(s,:)),real(u_hat(s,:)));
rnorm(real(gamma_mean_trace_smoother_sample(s,:)),real(u_hat(s,:)));



% reconstruction
Dim_Grid = 25;
ns=2*pi;

[xx,yy] = meshgrid(linspace(0,ns,Dim_Grid), linspace(0,ns,Dim_Grid));
x_vec = [reshape(xx,[],1), reshape(yy,[],1)]; 

return
figure() 
u0=[];
% for i = 2:10:round(T/dt/100)
for i=1:10:N
% for i=15000
i1=1+50*(i-1);
    u = exp(1i * x_vec * kk*(2*pi)/ns) * (u_hat(:,i) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk*(2*pi)/ns) * (u_hat(:,i) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    imagesc([0,ns],[0,ns],sqrt(real(u).^2+real(v).^2));colorbar
    hold on
    quiver(xx, yy, u, v, 'linewidth',1)
% hold on;
% 
    xlim([0, ns ])
    ylim([0, ns ])    
%     xlim([-pi, pi ])
%     ylim([-pi, pi ])
    box on    
    title(['t = ', num2str(i*dt)])
%     pause(.001);
hold on
scatter(xexact(1,i),yexact(1,i),'r')
pause()
% close all
    u0=[u0,sqrt(u(:).^2+v(:).^2)];
end
% max(u0(:))


















return
figure
for i = 1:4
    subplot(2,2,i)
    if i == 1
        plot(dt:dt:N*dt, real(u_hat(1,:)), 'b', 'linewidth',2)
        title(['(a) gravity mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
    elseif i == 2
        plot(dt:dt:N*dt, real(u_hat(Dim_Ug*2+1,:)), 'b', 'linewidth',2)
        title(['(b) GB mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
    elseif i == 3
        plot(dt:dt:N*dt, real(u_hat(end-1,:)), 'b', 'linewidth',2)
        title('(c) Zonal background flow','fontsize',14)
    elseif i == 4
        plot(dt:dt:N*dt, real(u_hat(end,:)), 'b', 'linewidth',2)
        title('(d) Meridional background flow','fontsize',14)
    end
    set(gca,'fontsize',12)
    box on
    xlabel('t')
end