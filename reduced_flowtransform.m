%%%% reduced impressible flow 
%%%% transform without localization but with variable number of tracers
clear;close all;rng(1);
dt=1/200;T=50;N=T/dt;regularizeval=10^(-3);gamma=1.5;
fprintf('gamma is %2.2f\n',gamma);
K_max=2;Kuse=2;L=12*3;diagcov=1;f0=0;if K_max~=Kuse;reduced=1;else reduced=0;end;
if Kuse>K_max;Kuse=K_max;reduced=0;end
E_0=1;k_0=2;alpha=3;3;nosmoother=1;backsample=0;
dimuhat=(2*K_max+1)^2;Kmax=K_max;
[ky,kx]=meshgrid([0:K_max,-K_max:-1],[0:K_max,-K_max:-1]);
kk=[kx(:),ky(:)]';
rk=zeros(size(kk));
for i=2:dimuhat
        ki=kk(:,i);
%         rk(:,i)=1i*[-ki(2);ki(1)]/norm(ki);
        rk(:,i)=1i*[-ki(2);ki(1)]/sqrt(ki(1)^2+ki(2)^2+1);
end
% ix=min(K_max,k_0);iy=min(K_max,k_0);ind=getind(ix,iy,kk);%% sample plot
if reduced
ix=k_0;iy=k_0;ind=getind(ix,iy,kk);%% sample plot    
else
ix=K_max;iy=K_max;ind=getind(ix,iy,kk);%% sample plot
end
redind=getreducedind(kk,K_max,Kuse);
fprintf('L is %d\n',L);
if reduced==0;redind=(1:dimuhat)';end;
rkred=rk(:,redind);kkred=kk(:,redind);
L1=2*pi;%%% size of the domain

sig_ex=0.25;
u_hat=zeros(dimuhat,N);
% % dkexact=-0.5*ones(dimuhat,1);
% % omegaexact0=0.5;omegaexact=formomega(dimuhat,omegaexact0,kk);
% % Fuexact=0.2*ones(dimuhat,1);
% % sigmauexact=0.5*ones(dimuhat,1);
% % sigmauhat=form_sigmatrix1(dimuhat,sigmauexact,kk);

[dkexact,omegaexact,Fuexact,sigmauhat,ene,sigmak2exact]=form_coeff(dimuhat,kk,E_0,k_0,alpha,Kmax,f0);
dkexact=-dkexact;

% % % sigmak2exact=sigmak2exact/sqrt(2);
fixedcov=diag((sigmak2exact.^2)./(-dkexact+sqrt(dkexact.^2+L/sig_ex/sig_ex.*sigmak2exact.^2)));
fixedcov(1,1)=1;fixedcov=fixedcov(redind,redind);
%% generate true signal
% % utest=zeros(size(u_hat));utest1=zeros(size(u_hat));
disp('generating true signal...')
for i=2:N
    u_hat(:,i)=u_hat(:,i-1)+(dkexact+1i*omegaexact).*u_hat(:,i-1)*dt+...
        Fuexact*dt+sqrt(dt)*sigmauhat*randn(dimuhat,1);
    if max(abs(real(u_hat(:,i))))>10^8
        disp('error, data blow up')
    end
    tempcode=u_hat(:,i);
% %     utest1(:,N-i+1)=tempcode;
% %     utest(:,i)=tempcode;
    
end
xexact=zeros(L,N);
yexact=xexact;
xexact(:,1)=L1*rand(L,1);
yexact(:,1)=L1*rand(L,1);

% return
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
Dim_Y=dimuhat;Dim_Yr=length(redind);
Dim_X=2*L;
gamma_mean0 = zeros(Dim_Yr,1);
gamma_cov0 = eye(Dim_Yr)*0.01; 
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
a0=Fuexact(redind);
a1=diag(dkexact(redind)+1i*omegaexact(redind));
b1=sigmauhat(redind,redind);
x=xexact;y=yexact;
invBoB = 1 /sig_ex/ sig_ex * eye(2*L); % inverse of the square of the observational noise
invBoBdiag=1 /sig_ex/ sig_ex * ones(2*L,1);
b1b1t=b1 * b1';
time=1;
I=eye(dimuhat);
allnobs=zeros(dimuhat,N);
 halfindgb=[];halfindgb1=[];
for iy=-Kmax:Kmax;
    for ix=1:Kmax;
%         ix=3;iy=2;
        i1=getind(ix,iy,kk);
        halfindgb=[halfindgb;i1(1)];
                i1=getind(-ix,-iy,kk);
        halfindgb1=[halfindgb1;i1(1)];
    end
end

for iy=1:Kmax;
    for ix=0
%         ix=3;iy=2;
        i1=getind(ix,iy,kk);
        halfindgb=[halfindgb;i1(1)];
                        i1=getind(-ix,-iy,kk);
        halfindgb1=[halfindgb1;i1(1)];
    end
end 
dimuhat0=dimuhat;
kk0=kk;xy=kk0'/(2*K_max+1)*2*pi;
rketa=[1./sqrt(kk0(1,:).^2 + kk0(2,:).^2+1),...
    1./sqrt(kk0(1,:).^2 + kk0(2,:).^2)/sqrt(2)./sqrt(kk0(1,:).^2 + kk0(2,:).^2 + 1).*(kk0(1,:).^2 + kk0(2,:).^2),...
    1./sqrt(kk0(1,:).^2 + kk0(2,:).^2)/sqrt(2)./sqrt(kk0(1,:).^2 + kk0(2,:).^2 + 1).*(kk0(1,:).^2 + kk0(2,:).^2)];
rketa(1+dimuhat0)=0;rketa(1+dimuhat0*2)=0;

invG=exp(1i * xy*kk0);%%% inverse shallow water matrix
G=exp(-1i  * xy*kk0)/(2*K_max+1)^2;%%%%% shallow water matrix from physical domain to F domain
% invfft2mat=invG;
invG1=invG .* (ones(length(xy),1) * rk(1,1:dimuhat0)); 
invG2=invG .* (ones(length(xy),1) * rk(2,1:dimuhat0)); 
invG3=invG .*(ones(length(xy),1) * rketa(1:dimuhat0) );

Gm=[G,G,G];
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
    
    A1=[G1;G2];A1=A1(:,redind);
    
% idx=[1:20,1+L:20+L];


  idmeshdrift=get_local_vdID((2*K_max+1),gamma,xexact(:,i),yexact(:,i),L1,0);  
% gamma_mean=gamma_mean0;gamma_cov=zeros(dimuhat,dimuhat);gamma_cov(1,1)=1;

idx=[idmeshdrift{1},idmeshdrift{1}+L];allnobs1(i)=length(idx)/2;
testcase=1;
if testcase==1
covph=zeros(3*dimuhat,3*dimuhat);
Lu=invG1;
    u = Lu*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lu*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
    covph(1:dimuhat,1:dimuhat) = Lu*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lu' - ((Lu*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lu*gamma_cov0 * A1(idx,:)')') * dt;  

Lu=invG2;
    v = Lu*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lu*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
     covph(1+dimuhat:2*dimuhat,1+dimuhat:2*dimuhat) = Lu*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lu' - ((Lu*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lu*gamma_cov0 * A1(idx,:)')') * dt;  


Lu=invG3;
    eta = Lu*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lu*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
    covph(1+2*dimuhat:3*dimuhat,1+2*dimuhat:3*dimuhat) = Lu*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lu' - ((Lu*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lu*gamma_cov0 * A1(idx,:)')') * dt;  
% % else %%% wrong
% %     Lu=[invG1;invG2;invG3];
% % temp = Lu*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
% %         (Lu*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
% % covph = Lu*(gamma_cov0 + (a1 * gamma_cov0 + ...
% % gamma_cov0 * a1' + b1b1t)*dt)*Lu' - ((Lu*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lu*gamma_cov0 * A1(idx,:)')') * dt;     
% % u=temp(1:dimuhat0);v=temp(1+dimuhat0:2*dimuhat0);eta=temp(1+2*dimuhat0:3*dimuhat0);
    
end
Gu=G*u;
Gv=G*v;
gamma_mean=gamma_mean0;
for id=2:dimuhat0
%     if length(idmeshdrift{id})==0%%%% of no observation then no update
% % %     disp('no obs')
%     else
    ka=kk0(1,id);kb=kk0(2,id);
    if rk(1,id)==0
      gamma_mean(id)=Gv(id)/rk(2,id);  
    else
       gamma_mean(id)=Gu(id)/rk(1,id);   
    end
%     end
end
gamma_mean(1)=0;

gamma_cov=Gm*(covph)*Gm';%%%% transform back covariance matrix
  
%     gamma_cov =fixedcov;gamma_cov0=gamma_cov;

    % save the posterior statistics
    gamma_mean_trace(redind,i) = gamma_mean;

    % update
    gamma_mean0 = gamma_mean;     
%     if diagcov
%     gamma_cov0 = diag(diag(gamma_cov));
%     else
    gamma_cov0 = gamma_cov;
%     end

     if nosmoother~=1%%%% record for smoother
         if diagcov
  gamma_cov_trace(:,i) = diag(gamma_cov0);             
         else
    gamma_cov_trace(:,:,i) = gamma_cov0;
         end
     end   
    if mod(i,1/dt)==0
    u = (Ga*  gamma_mean_trace(:,i)); 
    v = (Gb*  gamma_mean_trace(:,i)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity,  record filter'  )
      max(abs(imag(u)))
    end
    uexactfilter(:,time)=real(u);
    vexactfilter(:,time)=real(v);
     
        time=time+1;  
    end
end
fprintf('mean obs is %2.2f\n',mean(allnobs1(2:end)));
figure();plot(dt:dt:T,real(u_hat(ind,:)),'r','linewidth',1.5);hold on;
plot(dt:dt:T,real(gamma_mean_trace(ind,:)),'b','linewidth',1.5);
setgca(16);title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
xlabel('t','fontsize',16)
% rnorm(real(gamma_mean_trace(ind,:)),real(u_hat(ind,:)));
rmscc(real(gamma_mean_trace(ind(1),:)),real(u_hat(ind(1),:)),1);

% imagesc(reshape(uexact(:,5).^2+vexact(:,5).^2,40,40))
exact_gamma_mean_trace=gamma_mean_trace;
if nosmoother
return
end
exact_gamma_cov_trace=gamma_cov_trace;

