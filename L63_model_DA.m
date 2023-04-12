clear;close all;
rng(10);Linewid=2;fonsize=18;
% Generating a true signal of the SPEKF-M model
T = 500;
dt = 5e-3;
N = round(T/dt);
x_truth = zeros(1,N);
y_truth = zeros(1,N);
z_truth = zeros(1,N);
Tacf=15;range=floor(50/dt):floor(150/dt);
rho = 28;
sigma = 10;
beta = 8/3;

% sigma_x = 5;
% sigma_y = 5;
% sigma_z = 5;

sigma_x = 10;
sigma_y = 1;
sigma_z = 1;

for i = 2:N
    x_truth(:,i) = x_truth(:,i-1) + sigma * (y_truth(:,i-1) - x_truth(:,i-1)) * dt + sigma_x * sqrt(dt) * randn; 
    y_truth(:,i) = y_truth(:,i-1) + (x_truth(:,i-1) * (rho - z_truth(:,i-1)) - y_truth(:,i-1)) * dt + sigma_y * sqrt(dt) * randn;
    z_truth(:,i) = z_truth(:,i-1) + (x_truth(:,i-1) * y_truth(:,i-1) - beta * z_truth(:,i-1)) * dt + sigma_z * sqrt(dt) * randn;
end
% F_u_est = F_u_truth;
% F_g_est = F_g_truth;
% d_u_est = d_u_truth;
% d_v_est = d_v_truth;
% d_g_est = d_g_truth;
% c_est = c_truth;
% sigma_u_est = sigma_u_truth;
% sigma_v_est = sigma_v_truth;
% sigma_g_est = sigma_g_truth;

rho_est = 28;
sigma_est = 10;
beta_est = 8/3;

sigma_x_est = 5;
sigma_y_est = 5;
sigma_z_est = 5;

gamma_mean_trace = zeros(2,N);
gamma_cov_trace = zeros(4,N);


a1 = zeros(2,2);
b1 = zeros(2,2);
% load DyadGM;

invBoB = 1/sigma_x_est/sigma_x_est;


% if mod(kk,1000) ==1
% disp(kk)
% end
gamma_mean0 = [0;0];
gamma_cov0 = zeros(2,2);
for i = 2:N

    u0 = x_truth(i-1);
    u = x_truth(i);
    a1 = [-1, -u0; u0, -beta_est];
        
    a0 = [rho_est * u0; 0];
    b1 = [sigma_y_est, 0;0, sigma_z_est];
    A0 = -sigma_est * u0;
    A1 = [sigma_est, 0];


    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + ...
    (gamma_cov0 * A1') * invBoB * (u-u0 - A0*dt-A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' ...
        - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     


%     [gamma_cov,LL,GG] = care(a1'*dt-eye(1)/2 , A1'*sqrt(invBoB*dt), b1 * b1'*dt + gamma_cov0);
%     gamma_mean = (eye(1) - a1*dt + (gamma_cov * A1') * invBoB * A1 * dt )\((gamma_cov * A1') * invBoB * (u-u0-A0*dt) +gamma_mean0 +a0*dt);

%     gamma_trace(i) = gamma;
    gamma_mean_trace(:,i) = gamma_mean;
    gamma_cov_trace(:,i) = reshape(gamma_cov,4,[]);

    gamma_mean0 = gamma_mean;
    gamma_cov0 = gamma_cov;

end

gamma_mean0_smoother = gamma_mean0;
gamma_cov0_smoother = gamma_cov0;
gamma_mean_trace_smoother = zeros(2,N);
gamma_mean_trace_smoother(:,end)=gamma_mean0;

% return
for i=N:-1:2
    R=gamma_cov_trace(:,i);R=reshape(R,2,2);invR=R\speye(2);
    u0 = x_truth(i);
    a1 = [-1, -u0; u0, -beta_est];
        
    a0 = [rho_est * u0; 0];
    b1 = [sigma_y_est, 0;0, sigma_z_est];

    
    gamma_mean_smoother=gamma_mean0_smoother+dt*(-a0-a1*gamma_mean0_smoother+(b1 * b1')*invR*(gamma_mean_trace(:,i)-gamma_mean0_smoother));
    gamma_cov_smoother=gamma_cov0_smoother-dt*( (a1+(b1*b1')*invR)*gamma_cov0_smoother+...
        gamma_cov0_smoother*(a1'+b1*b1'*invR)-b1*b1' );
        gamma_mean_trace_smoother(:,i-1) = gamma_mean_smoother;

    gamma_mean0_smoother = gamma_mean_smoother;
    gamma_cov0_smoother = gamma_cov_smoother;
end


% plot2line( gamma_mean_trace(1,range),gamma_mean_trace_smoother(1,range))
plot2line( y_truth(range),gamma_mean_trace_smoother(1,range))
plot2line( z_truth(range),gamma_mean_trace_smoother(2,range))
% plot2line( y_truth(range),gamma_mean_trace(1,range))
% plot2line( z_truth(range),gamma_mean_trace(2,range))
% return
%% smoother
gamma_mean0_smoother_sample = gamma_mean0;
gamma_mean_trace_smoother_sample = zeros(2,N);
gamma_mean_trace_smoother_sample(:,end)=gamma_mean0;
for i=N:-1:2
    R=gamma_cov_trace(:,i);R=reshape(R,2,2);invR=R\speye(2);
    u0 = x_truth(i);
    a1 = [-1, -u0; u0, -beta_est];
        
    a0 = [rho_est * u0; 0];
    b1 = [sigma_y_est, 0;0, sigma_z_est];

    
    gamma_mean_smoother_sample=gamma_mean0_smoother_sample+dt*(-a0-a1*gamma_mean0_smoother_sample+...
        (b1 * b1')*invR*(gamma_mean_trace(:,i)-gamma_mean0_smoother_sample))+sqrt(dt)*b1*randn(2,1);
gamma_mean_trace_smoother_sample(:,i-1) = gamma_mean_smoother_sample;

    gamma_mean0_smoother_sample = gamma_mean_smoother_sample;
end
plot2line( z_truth(range),gamma_mean_trace_smoother_sample(2,range))
% plot2line( y_truth(range),gamma_mean_trace_smoother_sample(1,range))
%  plot2line( y_truth(range),gamma_mean_trace(1,range))
figure();[~,r1]=ksdensity(x_truth);plot(r1,ksdensity(x_truth),'r','Linewidth',Linewid);
return
range=1:N;
plotpdf(y_truth,gamma_mean_trace(1,:),gamma_mean_trace_smoother(1,:),gamma_mean_trace_smoother_sample(1,:),range)
plotpdf(z_truth,gamma_mean_trace(2,:),gamma_mean_trace_smoother(2,:),gamma_mean_trace_smoother_sample(2,:),range)

plotacf(y_truth,gamma_mean_trace(1,:),gamma_mean_trace_smoother(1,:),gamma_mean_trace_smoother_sample(1,:),1:N,Tacf,dt)
plotacf(z_truth,gamma_mean_trace(2,:),gamma_mean_trace_smoother(2,:),gamma_mean_trace_smoother_sample(2,:),1:N,Tacf,dt)
% rmse(gamma_mean_trace(1,:),y_truth);
% rmse(gamma_mean_trace_smoother(1,:),y_truth);
% rmse(gamma_mean_trace_smoother_sample(1,:),y_truth);
% rmse(gamma_mean_trace(2,:),z_truth);
% rmse(gamma_mean_trace_smoother(2,:),z_truth);
% rmse(gamma_mean_trace_smoother_sample(2,:),z_truth);

rnorm(gamma_mean_trace(1,:),y_truth);
rnorm(gamma_mean_trace_smoother(1,:),y_truth);
rnorm(gamma_mean_trace_smoother_sample(1,:),y_truth);
rnorm(gamma_mean_trace(2,:),z_truth);
rnorm(gamma_mean_trace_smoother(2,:),z_truth);
rnorm(gamma_mean_trace_smoother_sample(2,:),z_truth);
return
figure
subplot(3,1,2)
hold on
% curve1 = gamma_mean_trace(1,:) + sqrt(gamma_cov_trace(1,:));
% curve2 = gamma_mean_trace(1,:) - sqrt(gamma_cov_trace(1,:));
% inBetween = [curve1, fliplr(curve2)];
% x2 = [dt:dt:T, fliplr(dt:dt:T)];
% p1 = fill(x2, inBetween,  [0.8,0.8,0.8],'linestyle','none');
p2 = plot(dt:dt:T,y_truth,'b');
p3 = plot(dt:dt:T,gamma_mean_trace(1,:),'r');
% legend([p2,p3,p1],'Truth','Pred','Confidence interval')
set(gca,'fontsize',12);
box on
% title('y','fontsize',16)
xlim([0,50])
subplot(3,1,3)
hold on
% curve1 = gamma_mean_trace(1,:) + sqrt(gamma_cov_trace(1,:));
% curve2 = gamma_mean_trace(1,:) - sqrt(gamma_cov_trace(1,:));
% inBetween = [curve1, fliplr(curve2)];
% x2 = [dt:dt:T, fliplr(dt:dt:T)];
% p1 = fill(x2, inBetween,  [0.8,0.8,0.8],'linestyle','none');
p2 = plot(dt:dt:T,z_truth,'b');
p3 = plot(dt:dt:T,gamma_mean_trace(2,:),'r');
legend('Truth','assimiated')
% legend([p2,p3,p1],'Truth','Pred','Confidence interval')
set(gca,'fontsize',12);
box on
% title('z','fontsize',16)
xlim([0,50])
subplot(3,1,1)
plot(dt:dt:T,x_truth,'b')
set(gca,'fontsize',12);
box on
% title('x','fontsize',16)
xlim([0,50])