%%%% layered model, conditional gaussian
%%%% with direction (lx,ly)
%%%% almost the same fit_layered 
clear;close all;
rng(10);
T=500;dt=1/1000;N=T/dt;L1=2*pi;rk=1;
L=20;Luse=min(L,20);Ltotal=L;fprintf('L is %d \n',Luse);
testcase1=1;%%%1:gauss 2:strong non gauss
if testcase1==1
    load dkfitc1;load omegakfitc1;load Ffitc1;load sigmafitc1
else
  load dkfitc2;load omegakfitc2;load Ffitc2;load sigmafitc2
end
b1b1t=sigmafit * sigmafit';
fprintf('test case %d\n',testcase1);
testcase=1;
if testcase==1
constant_cov=0;diagcov=0;
elseif testcase==2
constant_cov=0;diagcov=1;
elseif testcase==3
constant_cov=1;diagcov=1;
end
Kmax=6;p=2;%theta=pi/3;
fprintf('p is %2.2f\n',p);
fprintf('test case is %d\n',testcase);
fprintf('L is %d \n',Luse);
kk0=(1:Kmax)';lx=0.8;ly=sqrt(1-lx^2);fprintf('(lx, ly) are (%2.2f,%2.2f)\n',lx,ly);
dk=0.0125*ones(Kmax,1);du=0.0125;H1=1;H2=1/2;
beta=2;sigmak=zeros(Kmax,1);hk=zeros(Kmax,1);
x=zeros(L,N);sig_ex=.25;
y=x;gap=1;dt0=1*dt;N0=T/dt0;
fprintf('Kmax is %d\n',Kmax);
fprintf('sig ex is %2.2f\n',sig_ex);
x(:,1)=L1*rand(L,1);
y(:,1)=L1*rand(L,1);
if testcase==1;
    sigmau=1/2/sqrt(2);sigmak(1)=1/1/sqrt(2);sigmak(2)=1/2/sqrt(2);
else
    sigmau=1/1/sqrt(2);sigmak(1)=1/2/sqrt(2)/2;sigmak(2)=1/2/sqrt(2)/2;
end
hk(1)=H1/2-H1/2*1i;hk(2)=H2/2-H2/2*1i;

% a=real(psik);b=imag(psik);c=sum(a.^2+b.^2,2);c(1)/c(2)
for ik=3:Kmax
    if testcase==1
    sigmak(ik)=1/ik^p/sqrt(2);
    else
     sigmak(ik)=1/ik^p/sqrt(2)/2;   
    end
    theta=rand(1)*2*pi;
    hk(ik)=sin(theta)/4/(ik)^p-cos(theta)/4/(ik)^p*1i;
    
end
hkstar=real(hk)-1i*imag(hk);
psik=zeros(Kmax,N);u=zeros(1,N);
tempu=u(:,1);
temppsik=psik(:,1);

kk1dfull=zeros(2,2*Kmax);
for i=1:Kmax
   kk1dfull(1,2*i-1)=i*lx;
   kk1dfull(1,2*i)=-i*lx;
   kk1dfull(2,2*i-1)=i*ly;
   kk1dfull(2,2*i)=-i*ly;   
end

kk1dfullf=[[0;0],kk1dfull];
repk1d=repmat(kk1dfull(1,:),L,1)*1i;
kx=kk1dfullf(1,:);ky=kk1dfullf(2,:);
dX = 1i*kx*2*pi/L1;
dY = 1i*ky*2*pi/L1;
dX=dX(:);dY=dY(:);

%%% generate data
rng(10);
dim=length(dkfit);
ufree=zeros(dim, N);
uold=zeros(dim,1);
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end    
        x_loc = [x(:,i-1),y(:,i-1)];  
   u1=uold+(-dkfit+1i*omegakfit).*uold*dt+Ffit*dt+sqrt(dt)*sigmafit*randn(dim,1);
  
    G0=zeros(L,2*Kmax+1);G0(:,1)=1;
    G1 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (-ones(L,1) * transpose(dY));
    G2 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (ones(L,1) * transpose(dX));
    G1=G1+G0;G2=G2+G0;
    ux=G1*uold;uy=G2*uold;
    

if max(abs(imag(ux(:))))>10^(-6) || max(abs(imag(uy(:))))>10^(-6)
    disp('error, complex velicity')
end
ux=real(ux);uy=real(uy);
 uold=u1;
%% update position
    x(:,i) = x(:,i-1) + ux * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in x
    y(:,i) = y(:,i-1) + uy * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in y
    x(:,i) = mod(x(:,i),L1);
    y(:,i) = mod(y(:,i),L1); 
     ufree(:,i)=u1;
end
%% conditional gaussian
dimuhat=2*Kmax+1;redind=1:dimuhat;nosmoother=1;
Dim_Y=dimuhat;Dim_Yr=length(redind);
Dim_X=2*L;
gamma_mean0 = zeros(Dim_Yr,1);
gamma_cov0 = eye(Dim_Yr)*0.01; 
if constant_cov
    gamma_cov0=fixedcov;
end
% gamma_cov0 =fixedcov;
gamma_mean_trace=zeros(Dim_Y,N);%%% store all
gamma_mean_trace(redind,1) = gamma_mean0;
if nosmoother~=1%%%% only store reduced dof
    if diagcov
gamma_cov_trace=zeros(Dim_Yr,N);
gamma_cov_trace(:,1) = diag(gamma_cov0);        
    else
gamma_cov_trace=zeros(Dim_Yr,Dim_Yr,N);
gamma_cov_trace(:,:,1) = gamma_cov0;
    end
end

A0=zeros(Dim_X,1);
a0=Ffit(redind);
a1=diag(-dkfit(redind)+1i*omegakfit(redind));


invBoB = 1 /sig_ex/ sig_ex * eye(2*L); % inverse of the square of the observational noise
invBoBdiag=1 /sig_ex/ sig_ex * ones(2*L,1);
 G0=zeros(L,2*Kmax+1);G0(:,1)=1;
time=1;
disp('data assimilation (filter)......')
for i = 2:N
    if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
    end
    % observational operator 
    x_loc = [x(:,i-1),y(:,i-1)];
   
    G1 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (-ones(L,1) * transpose(dY));
    G2 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (ones(L,1) * transpose(dX));
    G1=G1+G0;G2=G2+G0;
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
    
    A1=[G1;G2];A1=A1(:,redind);
    
% idx=[1:20,1+L:20+L];

    % run the data assimilation for posterior mean and posterior covariance
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * (invBoBdiag.*  (diff_xy - A0*dt-A1 * gamma_mean0 * dt));
if constant_cov
    gamma_cov=fixedcov;
else
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1b1t - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;         
end

    

    % save the posterior statistics
    gamma_mean_trace(redind,i) = gamma_mean;

    % update
    gamma_mean0 = gamma_mean;     
    if diagcov
    gamma_cov0 = diag(diag(gamma_cov));
    else
    gamma_cov0 = gamma_cov;
    end                


% %          if diagcov
% %   gamma_cov_trace(:,i) = diag(gamma_cov0);             
% %          else
% %     gamma_cov_trace(:,:,i) = gamma_cov0;
% %          end
end
i=1;rmscc(real( gamma_mean_trace(i,:)),real(ufree(i,:)),1);
i=2;rmscc(real( gamma_mean_trace(i,:)),real(ufree(i,:)),1);
i=4;rmscc(real( gamma_mean_trace(i,:)),real(ufree(i,:)),1);