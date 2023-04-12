%%%% reduced shallow water, no transform test
%%%% transform without localization but with variable number of tracers

clear;close all;rng(1);
dt=1/200;T=20;N=floor(T/dt);regularizeval=10^(-3);
gamma=2.0;fprintf('gamma is %2.2f\n',gamma);
K_max=1;Kuse=1;L=12*3;diagcov=0;if K_max~=Kuse;reduced=1;else reduced=0;end;
if Kuse>K_max;Kuse=K_max;reduced=0;end
eps=1;E_0=1;E_1=0.12;k_0=2;alpha=3;5/3;nosmoother=1;backsample=0;
fg=0;fr=0;
dimuhat0=(2*K_max+1)^2;dimuhat=dimuhat0*3;Kmax=K_max;
[ky,kx]=meshgrid([0:K_max,-K_max:-1],[0:K_max,-K_max:-1]);
kk=[kx(:),ky(:)]';

rk1 = [1./sqrt(kk(1,:).^2 + kk(2,:).^2+1) .* (-1i * kk(2,:));
    1./sqrt(kk(1,:).^2 + kk(2,:).^2+1) .* (1i * kk(1,:))];
rk1(:,1) = [0;0]; 
rk2 = [1./sqrt(kk(1,:).^2 + kk(2,:).^2)/sqrt(2)./sqrt(kk(1,:).^2 + kk(2,:).^2 + 1) .* (1i * kk(2,:) + kk(1,:) .* sqrt(kk(1,:).^2 + kk(2,:).^2 + 1));
    1./sqrt(kk(1,:).^2 + kk(2,:).^2)/sqrt(2)./sqrt(kk(1,:).^2 + kk(2,:).^2 + 1) .* (-1i * kk(1,:) + kk(2,:) .* sqrt(kk(1,:).^2 + kk(2,:).^2 + 1))];
% rk2(:,1) =[1i;1]/sqrt(2);
rk2(:,1) =[0;0]/sqrt(2);
rk3 = -[1./sqrt(kk(1,:).^2 + kk(2,:).^2)/sqrt(2)./sqrt(kk(1,:).^2 + kk(2,:).^2 + 1) .* (1i * kk(2,:) - kk(1,:) .* sqrt(kk(1,:).^2 + kk(2,:).^2 + 1));
    1./sqrt(kk(1,:).^2 + kk(2,:).^2)/sqrt(2)./sqrt(kk(1,:).^2 + kk(2,:).^2 + 1) .* (-1i * kk(1,:) - kk(2,:) .* sqrt(kk(1,:).^2 + kk(2,:).^2 + 1))];
% rk3(:,1) = [-1i;1]/sqrt(2);
rk3(:,1) = [0;0]/sqrt(2);

redind=getreducedind(kk,K_max,Kuse);redind=[redind;redind+dimuhat0;redind+2*dimuhat0];
kk0=kk;kk=[kk,kk,kk];negindkk=zeros(dimuhat0,1);
for i=1:length(kk0)
ix=kk0(1,i);iy=kk0(2,i);[a]=find(kk0(1,:)==-ix);[b]=find(kk0(2,:)==-iy);negindkk(i)=intersect(a,b); 
end
% ix=min(K_max,k_0);iy=min(K_max,k_0);ind=getind(ix,iy,kk);%% sample plot
ix=K_max;iy=K_max;ind=getind(ix,iy,kk);%% sample plot
rk=zeros(2,dimuhat);rk(:,1:2*dimuhat0)=[rk1,rk2];
rk(:,negindkk+2*dimuhat0)=rk3;%%% guarantee the symmetry
if reduced==0;redind=(1:dimuhat)';end;
rkred=rk(:,redind);kkred=kk(:,redind);
L1=2*pi;%%% size of the domain
xexact=zeros(L,N);
yexact=xexact;

sig_ex=0.5;
u_hat=zeros(dimuhat,N);
% % dkexact=-0.5*ones(dimuhat,1);
% % omegaexact0=0.5;omegaexact=formomega(dimuhat,omegaexact0,kk);
% % Fuexact=0.2*ones(dimuhat,1);
% % sigmauexact=0.5*ones(dimuhat,1);
% % sigmauhat=form_sigmatrix1(dimuhat,sigmauexact,kk);

[dkexactb,omegaexactb,Fuexactb,sigmauhatb,eneb,sigmak2exactb]=form_coeff(dimuhat0,kk0,E_0,k_0,alpha,Kmax,fg);%%%%gb
[dkexactg,omegaexactg,Fuexactg,sigmauhatg,eneg,sigmak2exactg]=form_coeff(dimuhat0,kk0,E_1,k_0,alpha,Kmax,fr);%%%%gravity
dkexact=-[dkexactb;dkexactg;dkexactg];
omegaexact=[zeros(dimuhat0,1);1/eps*sqrt(kk0(1,:).^2 + kk0(2,:).^2)';-1/eps*sqrt(kk0(1,negindkk).^2 + kk0(2,negindkk).^2)'];
Fuexact=[Fuexactb;Fuexactg;Fuexactg];
sigmauhat=zeros(dimuhat,dimuhat);sigmauhat(1:dimuhat0,1:dimuhat0)=sigmauhatb;
sigmauhat(1+dimuhat0:end,1+dimuhat0:end)=form_gravity_sigma(kk0,sigmak2exactg);
% return
% % % sigmak2exact=sigmak2exact/sqrt(2);
% fixedcov=diag((sigmak2exact.^2)./(-dkexact+sqrt(dkexact.^2+L/sig_ex/sig_ex.*sigmak2exact.^2)));
% fixedcov(1,1)=1;fixedcov=fixedcov(redind,redind);
%% generate true signal
% % utest=zeros(size(u_hat));utest1=zeros(size(u_hat));
disp('generating true signal...')
for i=2:N
    if mod(i,floor(N/5)) == 0;fprintf('rest step is %d\n',N-i);end
    u_hat(:,i)=u_hat(:,i-1)+(dkexact+1i*omegaexact).*u_hat(:,i-1)*dt+...
        Fuexact*dt+1*sqrt(dt)*sigmauhat*randn(dimuhat,1);
    if max(abs(real(u_hat(:,i))))>10^8
        disp('error, data blow up')
    end
    tempcode=u_hat(:,i);
% %     utest1(:,N-i+1)=tempcode;
% %     utest(:,i)=tempcode;
    
end
xexact(:,1)=L1*rand(L,1);
yexact(:,1)=L1*rand(L,1);
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
disp('generating obsveration...')
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

xy=kk0'/(2*K_max+1)*2*pi;
rketa=[1./sqrt(kk0(1,:).^2 + kk0(2,:).^2+1),...
    1./sqrt(kk0(1,:).^2 + kk0(2,:).^2)/sqrt(2)./sqrt(kk0(1,:).^2 + kk0(2,:).^2 + 1).*(kk0(1,:).^2 + kk0(2,:).^2),...
    1./sqrt(kk0(1,:).^2 + kk0(2,:).^2)/sqrt(2)./sqrt(kk0(1,:).^2 + kk0(2,:).^2 + 1).*(kk0(1,:).^2 + kk0(2,:).^2)];
rketa(1+dimuhat0)=0;rketa(1+dimuhat0*2)=0;
    invG1 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rk(1,:))); % Fourier bases for u
    invG2 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rk(2,:))); % Fourier bases for v
    invG3 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rketa)); % Fourier bases for v
%     u = (invG1* data); 
%     v = (invG2* data); 
%     eta=(invG3* data);
    G=exp(-1i  * xy*kk0)/(2*K_max+1)^2;%%%%% shallow water matrix from physical domain to F domain
    Gm=[G,zeros(dimuhat0,dimuhat0),zeros(dimuhat0,dimuhat0);
        zeros(dimuhat0,dimuhat0),G,zeros(dimuhat0,dimuhat0);
        zeros(dimuhat0,dimuhat0),zeros(dimuhat0,dimuhat0),G];

    
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
    


    % run the data assimilation for posterior mean and posterior covariance
%     gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (diff_xy - A0*dt-A1 * gamma_mean0 * dt);
%     gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     

if i==2
  idmeshdrift=get_local_vdID((2*K_max+1),gamma,xexact(:,i),yexact(:,i),L1,0);  
  i
end 

idx=[idmeshdrift{1},idmeshdrift{1}+L];
% idx=1:2*L;
% idx=[1:12,(1:12)+L];
allnobs1(i)=length(idx)/2;
testcase=2;
if testcase==1
    
covph=zeros(3*dimuhat0,3*dimuhat0);
Lu=invG1;
u = Lu*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lu*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
covph(1:dimuhat0,1:dimuhat0) = Lu*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lu' - ((Lu*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lu*gamma_cov0 * A1(idx,:)')') * dt;     
Lu=invG2;
v = Lu*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lu*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
covph(1+dimuhat0:2*dimuhat0,1+dimuhat0:2*dimuhat0) = Lu*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lu' - ((Lu*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lu*gamma_cov0 * A1(idx,:)')') * dt;     
Lu=invG3;
eta = Lu*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lu*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
covph(1+2*dimuhat0:3*dimuhat0,1+2*dimuhat0:3*dimuhat0) = Lu*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lu' - ((Lu*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lu*gamma_cov0 * A1(idx,:)')') * dt;     
else
Lu=[invG1;invG2;invG3];
temp = Lu*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lu*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
covph = Lu*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lu' - ((Lu*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lu*gamma_cov0 * A1(idx,:)')') * dt;     
u=temp(1:dimuhat0);v=temp(1+dimuhat0:2*dimuhat0);eta=temp(1+2*dimuhat0:3*dimuhat0);
end
  ub=G*u;
  vb=G*v;
  etab=G*eta;
  
  gamma_mean=gamma_mean0;
    for ii=2:length(kk0)
       ind=[ii,ii+dimuhat0,ii+2*dimuhat0];
       ux=ub(ii);vx=vb(ii);etax=etab(ii);
       rmatrix=[rk(:,ind);rketa(:,ind)];
       gamma_mean(ind)=rmatrix\[ux;vx;etax];
    end 
gamma_cov0temp=Gm*(covph)*Gm';
% gamma_cov=gamma_cov0temp;
gamma_cov=zeros(size(gamma_cov0temp));
    for ii=2:length(kk0)
       ind=[ii,ii+dimuhat0,ii+2*dimuhat0];
       ux=ub(ii);vx=vb(ii);etax=etab(ii);
       rmatrix=[rk(:,ind);rketa(:,ind)];
      gamma_cov(ind,ind)=inv(rmatrix)*gamma_cov0temp(ind,ind)*inv(rmatrix');                  
    end
    gamma_cov([1,1+dimuhat0,1+2*dimuhat0],[1,1+dimuhat0,1+2*dimuhat0])=eye(3);

    
%     gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * (invBoBdiag.*  (diff_xy - A0*dt-A1 * gamma_mean0 * dt));
%     gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1b1t - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     

    
    % gamma_cov =fixedcov;gamma_cov0=gamma_cov;

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
        gamma_cov_trace(:,:,i) = gamma_cov0;
    if mod(i,1/dt)==0
    u = (Ga*  gamma_mean_trace(:,i)); 
    v = (Gb*  gamma_mean_trace(:,i)); 
    uexactfilter(:,time)=real(u);
    vexactfilter(:,time)=real(v);
     
        time=time+1;  
    end
    
     if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity,  record filter'  )
%       i
      max(abs(imag(u)))
    end
end
fprintf('mean obs is %2.2f\n',mean(allnobs1(2:end)));

ind0=ind(2);gap=200;Dt=gap*dt;
figure();subplot(2,1,1)
plot(dt:dt*gap:T,real(u_hat(ind(2),1:gap:end)),'r','linewidth',1.5);hold on;
plot(dt:dt*gap:T,real(gamma_mean_trace(ind(2),1:gap:end)),'b','linewidth',1.5);
setgca(16);title(['Filter, Gravity mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
subplot(2,1,2)
plot(dt:dt*gap:T,real(u_hat(ind(1),1:gap:end)),'r','linewidth',1.5);hold on;
plot(dt:dt*gap:T,real(gamma_mean_trace(ind(1),1:gap:end)),'b','linewidth',1.5);
setgca(16);title(['Filter, GB mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
xlabel('t','fontsize',16)

% rnorm(real(gamma_mean_trace(ind(1),:)),real(u_hat(ind(1),:)));
% rnorm(real(gamma_mean_trace(ind(2),:)),real(u_hat(ind(2),:)));
rmscc(real(gamma_mean_trace(ind(1),:)),real(u_hat(ind(1),:)),1);
rmscc(real(gamma_mean_trace(ind(2),:)),real(u_hat(ind(2),:)),1);
% return
% imagesc(reshape(uexact(:,5).^2+vexact(:,5).^2,40,40))
exact_gamma_mean_trace=gamma_mean_trace;

