clear;
close all
MM=2000;
dt = 5e-3;
Linewid=2;fonsize=18;
T = 200;Ti=5;Ts=10;
Tuv1=zeros(Ti/dt,T-Ts);Tvu1=zeros(Ti/dt,T-Ts);%%%% exact
Tuv2=zeros(Ti/dt,T-Ts);Tvu2=zeros(Ti/dt,T-Ts);%%% filter
Tuv3=zeros(Ti/dt,T-Ts);Tvu3=zeros(Ti/dt,T-Ts);%%%%% smoother

N = round(T/dt);
u=zeros(N,1);
v=zeros(N,1);
range=floor(155/dt):floor(175/dt);
Tacf=5;
F_u=1;d_v=0.8;d_u=0.8;sigma_v=2;c=1.2;sigma_u=0.5;

for i=2:N
    u(i)=u(i-1)+dt*((-d_u+c*v(i-1))*u(i-1)+F_u)+sigma_u*sqrt(dt)*randn;
    v(i)=v(i-1)+(-d_v*v(i-1)-c*u(i-1)^2)*dt+sigma_v*sqrt(dt)*randn;
end
% figure();plot(u(range))
% figure();plot(v(range ))
% return
gamma_mean_trace=zeros(1,N);
gamma_cov_trace=zeros(1,N);
gamma_mean0=0;
gamma_cov0=0;
invBoB=1/sigma_u/sigma_u;
for i = 2:N

    u0 = u(i-1);
    ua = u(i);
    a0 = -c*u0^2;
        
    a1 = -d_v;
    b1 = sigma_v;
    A0 = -d_u*u0+F_u;
    A1 = c*u0;



    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + ...
    (gamma_cov0 * A1') * invBoB * (ua-u0 - A0*dt-A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' ...
        - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     


    gamma_cov_trace(i) = gamma_cov;
    gamma_mean_trace(:,i) = gamma_mean;

    gamma_mean0 = gamma_mean;
    gamma_cov0 = gamma_cov;

end
% figure();plot(v(floor(155/dt):floor(175/dt) ))
plot2line(v(floor(155/dt):floor(175/dt) ),gamma_mean_trace(floor(155/dt):floor(175/dt) ))
% return
gamma_mean0_smoother = gamma_mean0;
gamma_cov0_smoother = gamma_cov0;
gamma_mean_trace_smoother = zeros(1,N);
gamma_mean_trace_smoother (end)=gamma_mean;
% gamma_cov_trace_smoother = zeros(1,N);
% return
for i=N:-1:2
    R=gamma_cov_trace(:,i);
    
    u0 = u(i);
   a0=-c*u0^2;  
   a1 = -d_v;b1 = sigma_v;
    
    gamma_mean_smoother=gamma_mean0_smoother+dt*(-a0-a1*gamma_mean0_smoother+...
        (b1 * b1')*inv(R)*(gamma_mean_trace(i)-gamma_mean0_smoother));
    gamma_cov_smoother=gamma_cov0_smoother-dt*( (a1+(b1*b1')*inv(R))*gamma_cov0_smoother+...
        gamma_cov0_smoother*(a1'+b1*b1'*inv(R))-b1*b1' );
        gamma_mean_trace_smoother(:,i-1) = gamma_mean_smoother;


    gamma_mean0_smoother = gamma_mean_smoother;
    gamma_cov0_smoother = gamma_cov_smoother;
end
% plot2line( v(floor(155/dt):floor(175/dt)),gamma_mean_trace(floor(155/dt):floor(175/dt)))
plot2line( v(floor(155/dt):floor(175/dt)),gamma_mean_trace_smoother(floor(155/dt):floor(175/dt)))

parfor j=Ts:T
    
[uv,vu]=dyad_infolocal(u(j/dt),v(j/dt), Ti, MM,c,F_u,sigma_u,d_v,sigma_v,dt,d_u);    
 Tuv1(:,j)=uv;   Tvu1(:,j)=vu;
[uv,vu]=dyad_infolocal(u(j/dt),gamma_mean_trace(j/dt), Ti, MM,c,F_u,sigma_u,d_v,sigma_v,dt,d_u);     
Tuv2(:,j)=uv;   Tvu2(:,j)=vu;
[uv,vu]=dyad_infolocal(u(j/dt),gamma_mean_trace_smoother(j/dt), Ti, MM,c,F_u,sigma_u,d_v,sigma_v,dt,d_u);     
Tuv3(:,j)=uv;   Tvu3(:,j)=vu;
end
Tuv1(1,:)=[];Tuv2(1,:)=[];Tuv3(1,:)=[];
Tvu1(1,:)=[];Tvu2(1,:)=[];Tvu3(1,:)=[];
Tuv1(:,1:Ts-1)=[];Tuv2(:,1:Ts-1)=[];Tuv3(:,1:Ts-1)=[];
Tvu1(:,1:Ts-1)=[];Tvu2(:,1:Ts-1)=[];Tvu3(:,1:Ts-1)=[];
return
figure();
subplot(1,3,1);
imagesc(Tuv1);set(gca,'ydir','normal');colorbar;title('exact, u to v');
subplot(1,3,2);
imagesc(Tuv2);set(gca,'ydir','normal');colorbar;title('filter, u to v');
subplot(1,3,3);
imagesc(Tuv3);set(gca,'ydir','normal');colorbar;title('smoother, u to v');


figure();
subplot(1,3,1);
imagesc(Tvu1);set(gca,'ydir','normal');colorbar;title('exact, v to u');
subplot(1,3,2);
imagesc(Tvu2);set(gca,'ydir','normal');colorbar;title('filter, v to u');
subplot(1,3,3);
imagesc(Tvu3);set(gca,'ydir','normal');colorbar;title('smoother, v to u');

save Tuv1 Tuv1
save Tuv2 Tuv2
save Tuv3 Tuv3
save Tvu1 Tvu1
save Tvu2 Tvu2
save Tvu3 Tvu3

load Tuv1 
load Tuv2 
load Tuv3  
load Tvu1 
load Tvu2 
load Tvu3
% get Tuv1.mat Tuv2.mat Tuv3.mat Tvu1.mat Tvu2.mat Tvu3.mat
return
%% backward
gamma_mean0_smoother_sample = gamma_mean;
gamma_mean_trace_smoother_sample = zeros(1,N);
for i=N:-1:2
    R=gamma_cov_trace(:,i);
    u0 = u(i);

    a1 = -d_v;
    b1 = sigma_v;

     a0 = -c*u0^2;
   gamma_mean_smoother_back=gamma_mean0_smoother_sample+...
       dt*(-a0-a1*gamma_mean0_smoother_sample)+dt*(b1 * b1')*inv(R)*(gamma_mean_trace(i)-gamma_mean0_smoother_sample)+b1*sqrt(dt)*randn;
 gamma_mean_trace_smoother_sample(:,i-1) = gamma_mean_smoother_back;
 gamma_mean0_smoother_sample = gamma_mean_smoother_back;
end

% range=1:N;
plot2line( v(range),gamma_mean_trace_smoother_sample(range))


