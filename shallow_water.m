%%%% reduced shallow water
clear;close all;rng(1);
dt=1/500;T=100;N=floor(T/dt);regularizeval=10^(-3);
K_max=3;Kuse=K_max;if K_max~=Kuse;reduced=1;else reduced=0;end;fprintf('Kmax is K_max %d\n',K_max);
L=160;Luse=min(L,12);
testcase=2;nosmoother=1;
if testcase==1
constant_cov=0;diagcov=0;
elseif testcase==2
constant_cov=0;diagcov=1;
else
constant_cov=1;diagcov=1;
end
fprintf('test case is %\n',testcase);
if constant_cov==1;diagcov=1;end;if diagcov==0;constant_cov=0;end;
fprintf('constant_cov and diagcov are %d %d\n',constant_cov,diagcov);
if Kuse>K_max;Kuse=K_max;reduced=0;end
gap=1/dt;Dt=gap*dt;Np=T/Dt;

fprintf('L is %d\n',Luse);
eps=1;E_0=1;E_1=0.25;k_0=2;alpha=3;5/3;backsample=0;
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
sigmak2exact=[sigmak2exactb;sigmak2exactg;sigmak2exactg];
fixedcov=diag((sigmak2exact.^2)./(-dkexact+sqrt(dkexact.^2+Luse/sig_ex/sig_ex.*sigmak2exact.^2)));
fixedcov(1,1)=1;fixedcov(1+dimuhat0,1+dimuhat0)=1;
fixedcov(1+dimuhat0*2,1+2*dimuhat0)=1;fixedcov=fixedcov(redind,redind);

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
Dim_Grid = 30;
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

end
Rexact=zeros(dimuhat,dimuhat);
for j=1:dimuhat
% streamf=u_hat(j,:);
% Rexact(j,j)=var(streamf,1);
Rexact(j,j)=-sigmak2exact(j)/2/dkexact(j);
end
Rexact(1,1)=1;Rexact(1+dimuhat0,1+dimuhat0)=1;Rexact(1+dimuhat0*2,1+dimuhat0*2)=1;
xexact=xexact(1:Luse,:);
yexact=yexact(1:Luse,:);
L=Luse;
% return
%% filter
Dim_Y=dimuhat;Dim_Yr=length(redind);
Dim_X=2*L;
gamma_mean0 = zeros(Dim_Yr,1);
% gamma_cov0 = eye(Dim_Yr)*0.01; 
gamma_cov0 =fixedcov;
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

    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * (invBoBdiag.*  (diff_xy - A0*dt-A1 * gamma_mean0 * dt));
 if constant_cov
    gamma_cov=fixedcov;
else  
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1b1t - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
 end
    
    % gamma_cov =fixedcov;gamma_cov0=gamma_cov;

    % save the posterior statistics
    gamma_mean_trace(redind,i) = gamma_mean;

    % update
    gamma_mean0 = gamma_mean;     
         if nosmoother%%% only filter, test 3 cases
    if diagcov
    gamma_cov0 = diag(diag(gamma_cov));
    else
    gamma_cov0 = gamma_cov;
    end                
         else
         
        gamma_cov0 = gamma_cov; 
  
        end

%      if nosmoother~=1%%%% record for smoother
         if diagcov
  gamma_cov_trace(:,i) = diag(gamma_cov0);             
         else
    gamma_cov_trace(:,:,i) = gamma_cov0;
         end
%      end   
% %         gamma_cov_trace(:,:,i) = gamma_cov0;


    
end

ind0=ind(2);
% figure();subplot(2,1,1)
% plot(dt:dt*gap:T,real(u_hat(ind(2),1:gap:end)),'r','linewidth',1.5);hold on;
% plot(dt:dt*gap:T,real(gamma_mean_trace(ind(2),1:gap:end)),'b','linewidth',1.5);
% setgca(16);title(['Filter, Gravity mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
% subplot(2,1,2)
% plot(dt:dt*gap:T,real(u_hat(ind(1),1:gap:end)),'r','linewidth',1.5);hold on;
% plot(dt:dt*gap:T,real(gamma_mean_trace(ind(1),1:gap:end)),'b','linewidth',1.5);
% setgca(16);title(['Filter, GB mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
% xlabel('t','fontsize',16)

% rnorm(real(gamma_mean_trace(ind(1),:)),real(u_hat(ind(1),:)));
% rnorm(real(gamma_mean_trace(ind(2),:)),real(u_hat(ind(2),:)));
rmscc(real(gamma_mean_trace(ind(1),:)),real(u_hat(ind(1),:)),1);
rmscc(real(gamma_mean_trace(ind(2),:)),real(u_hat(ind(2),:)),1);
% return
% imagesc(reshape(uexact(:,5).^2+vexact(:,5).^2,40,40))
exact_gamma_mean_trace=gamma_mean_trace;
if nosmoother
 u_hatgb=u_hat;u_hatgb(dimuhat0+1:end,:)=0;
 filgb=exact_gamma_mean_trace; filgb(dimuhat0+1:end,:)=0;
  u_hatgr=u_hat;u_hatgr(1:dimuhat0,:)=0;
 filgr=exact_gamma_mean_trace; filgr(1:dimuhat0,:)=0;

 [rmstimefil,rmsfil,ccfil ]=compute_velerrorGBsw(u_hat(:,gap:gap:end),exact_gamma_mean_trace(:,gap:gap:end),rk,kk,L1,Dim_Grid,3);
 errfil=[mean(rmstimefil(3:end));mean(rmsfil);mean(ccfil)];
 
 fildata=[u_hat(:,gap:gap:end);exact_gamma_mean_trace(:,gap:gap:end)];

  [rmstimefilgb,rmsfilgb,ccfilgb ]=compute_velerrorGBsw(u_hatgb(:,gap:gap:end), filgb(:,gap:gap:end),rk,kk,L1,Dim_Grid,3);
 errfilgb=[mean(rmstimefilgb(3:end));mean(rmsfilgb);mean(ccfilgb)];
  [rmstimefilgr,rmsfilgr,ccfilgr ]=compute_velerrorGBsw(u_hatgr(:,gap:gap:end), filgr(:,gap:gap:end),rk,kk,L1,Dim_Grid,3);
 errfilgr=[mean(rmstimefilgr(3:end));mean(rmsfilgr);mean(ccfilgr)];
aa1=1:dimuhat0;aa2=1+dimuhat0:dimuhat0*3;
 if diagcov~=1
  [sigfilgb,disfilgb]=computeRentropy(u_hatgb(aa1,gap:gap:end),filgb(aa1,gap:gap:end),gamma_cov_trace(aa1,aa1,gap:gap:end),Rexact(aa1,aa1),diagcov);
  [sigfilgr,disfilgr]=computeRentropy(u_hatgr(aa2,gap:gap:end),filgr(aa2,gap:gap:end),gamma_cov_trace(aa2,aa2,gap:gap:end),Rexact(aa2,aa2),diagcov);
 [sigfil,disfil]=computeRentropy(u_hat(:,gap:gap:end),gamma_mean_trace(:,gap:gap:end),gamma_cov_trace(:,:,gap:gap:end),Rexact,diagcov);

 else
  [sigfilgb,disfilgb]=computeRentropy(u_hatgb(aa1,gap:gap:end),filgb(aa1,gap:gap:end),gamma_cov_trace(aa1,gap:gap:end),Rexact(aa1,aa1),diagcov);    
 [sigfilgr,disfilgr]=computeRentropy(u_hatgr(aa2,gap:gap:end),filgr(aa2,gap:gap:end),gamma_cov_trace(aa2,gap:gap:end),Rexact(aa2,aa2),diagcov);    
  [sigfil,disfil]=computeRentropy(u_hat(:,gap:gap:end),gamma_mean_trace(:,gap:gap:end),gamma_cov_trace(:,gap:gap:end),Rexact,diagcov);

 end
end
fprintf('L is %d\n',L);
if nosmoother
    
return
end
exact_gamma_cov_trace=gamma_cov_trace;
%% smoother
disp('smoother....')
gamma_mean0_smoother = gamma_mean0;
gamma_cov0_smoother = gamma_cov0;
gamma_mean_trace_smoother = zeros(Dim_Y,N);
% gamma_mean_trace_smoother(redind,end)=gamma_mean0;
gamma_mean_trace_smoother(redind,1)=gamma_mean0;%%% forward recording
gamma_mean_trace_red=gamma_mean_trace(redind,:);

time=1;
for i=N:-1:2
   if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',i);
   end
    j=N-i+2;%%% for forward recording
    % observational operator 
   if diagcov
        R=gamma_cov_trace(:,i);R=regularize(diag(R),regularizeval);
       invR=diag(1./diag(R));
   else
    R=gamma_cov_trace(:,:,i);R=regularize(R,regularizeval);
    invR=R\speye(Dim_Yr);
   end
% %     gamma_mean_smoother=gamma_mean0_smoother+dt*(-a0-a1*gamma_mean0_smoother+(b1 * b1')*invR*(gamma_mean_trace_red(:,i)-gamma_mean0_smoother));
% %     gamma_cov_smoother=gamma_cov0_smoother-dt*( (a1+(b1*b1')*invR)*gamma_cov0_smoother+...
% %         gamma_cov0_smoother*(a1'+b1*b1'*invR)-b1*b1' );

gamma_mean_smoother=gamma_mean0_smoother+dt*(-a0-a1*gamma_mean0_smoother+(b1b1t)*invR*(gamma_mean_trace_red(:,i)-gamma_mean0_smoother));    
% %     gamma_mean_trace_smoother(redind,i-1) = gamma_mean_smoother;%%%%expensive
gamma_mean_trace_smoother(redind,j) = gamma_mean_smoother;%%%%% forward record

    gamma_mean0_smoother = gamma_mean_smoother;
% %     gamma_cov0_smoother = gamma_cov_smoother;

% %     if mod(i,1/dt)==0
% %     u = (Ga*  gamma_mean_trace_smoother(:,i-1)); 
% %     v = (Gb*  gamma_mean_trace_smoother(:,i-1)); 
% %     if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
% %       disp('complex velocity, record smoother'  )
% %       max(abs(imag(u)))
% %     end
% %     uexactbackward(:,time)=real(u);
% %     vexactbackward(:,time)=real(v);
% %      
% %         time=time+1;  
% %     end      


    
end
gamma_mean_trace_smoother=gamma_mean_trace_smoother(:,end:-1:1);%%% reverse the record
% disp('smoother sw wave error');
%  [rmstimesm,rmsallsm,ccallsm ]=compute_velerrorGBsw(u_hat(:,1:gap:end),gamma_mean_trace_smoother(:,1:gap:end),rk,kk,L1,Dim_Grid,3);

ind0=ind(2);gap=200;Dt=gap*dt;
% figure();subplot(2,1,1)
% plot(dt:dt*gap:T,real(u_hat(ind(2),1:gap:end)),'r','linewidth',1.5);hold on;
% plot(dt:dt*gap:T,real(gamma_mean_trace_smoother(ind(2),1:gap:end)),'b','linewidth',1.5);
% setgca(16);title(['Smoother, Gravity mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
% subplot(2,1,2)
% plot(dt:dt*gap:T,real(u_hat(ind(1),1:gap:end)),'r','linewidth',1.5);hold on;
% plot(dt:dt*gap:T,real(gamma_mean_trace_smoother(ind(1),1:gap:end)),'b','linewidth',1.5);
% setgca(16);title(['Smoother, GB mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
% xlabel('t','fontsize',16)

rmscc(real(gamma_mean_trace_smoother(ind(1),:)),real(u_hat(ind(1),:)),1);
rmscc(real(gamma_mean_trace_smoother(ind(2),:)),real(u_hat(ind(2),:)),1);% rnorm(real(gamma_mean_trace(ind0,:)),real(u_hat(ind0,:)));
disp('...')
% rmscc(real(gamma_mean_trace(ind(1),:)),real(u_hat(ind(1),:)),1);
% rmscc(real(gamma_mean_trace(ind(2),:)),real(u_hat(ind(2),:)),1);
exact_gamma_mean_trace_smoother=gamma_mean_trace_smoother;
% disp('...')
disp('exact parameter smoother sw error');
 [rmstimeexactsm,rmsexactsm,ccexactsm  ]=compute_velerrorGBsw(u_hat(:,gap:gap:end),exact_gamma_mean_trace_smoother(:,gap:gap:end),rk,kk,L1,Dim_Grid,3);
 if diagcov~=1
 [sigsm,dissm]=computeRentropy(u_hat(:,gap:gap:end),exact_gamma_mean_trace_smoother(:,gap:gap:end),exact_gamma_cov_trace(:,:,gap:gap:end),Rexact,diagcov);
 else
 [sigsm,dissm]=computeRentropy(u_hat(:,gap:gap:end),exact_gamma_mean_trace_smoother(:,gap:gap:end),exact_gamma_cov_trace(:,gap:gap:end),Rexact,diagcov);    
 end
disp('***********************************************************************************************************************');
if backsample==0
    return
end
%% backward sampling
disp('backward sampling...')
time=1;
gamma_mean0_smoother_sample = gamma_mean0;
gamma_mean_trace_smoother_sample = zeros(Dim_Y,N);
% gamma_mean_trace_smoother_sample(redind,end)=gamma_mean0;
gamma_mean_trace_smoother_sample(redind,1)=gamma_mean0;
for i=N:-1:2
   if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',i);
   end   
    j=N-i+2;
   if diagcov
      R=gamma_cov_trace(:,i);R=regularize(diag(R),regularizeval);
       invR=diag(1./diag(R));
   else
    R=gamma_cov_trace(:,:,i);R=regularize(R,regularizeval);
    invR=R\speye(Dim_Yr);
   end

% %     gamma_mean_smoother_sample=gamma_mean0_smoother_sample+dt*(-a0-a1*gamma_mean0_smoother_sample+...
% %         (b1 * b1')*invR*(gamma_mean_trace_red(:,i)-gamma_mean0_smoother_sample))+sqrt(dt)*b1*randn(Dim_Yr,1);
% % gamma_mean_trace_smoother_sample(redind,i-1) = gamma_mean_smoother_sample;

    gamma_mean_smoother_sample=gamma_mean0_smoother_sample+dt*(-a0-a1*gamma_mean0_smoother_sample+...
        (b1b1t)*invR*(gamma_mean_trace_red(:,i)-gamma_mean0_smoother_sample))+sqrt(dt)*b1*randn(Dim_Yr,1);
gamma_mean_trace_smoother_sample(redind,j) = gamma_mean_smoother_sample;

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
gamma_mean_trace_smoother_sample=gamma_mean_trace_smoother_sample(:,end:-1:1);%%% reverse the record

figure();plot(dt:dt:T,real(u_hat(ind,:)),'r','linewidth',1.5);hold on;
plot(dt:dt:T,real(gamma_mean_trace_smoother_sample(ind,:)),'b','linewidth',1.5);
setgca(16);title(['Backward sampling, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
xlabel('t','fontsize',16)

rnorm(real(gamma_mean_trace(ind,:)),real(u_hat(ind,:)));
rnorm(real(gamma_mean_trace_smoother(ind,:)),real(u_hat(ind,:)));
rnorm(real(gamma_mean_trace_smoother_sample(ind,:)),real(u_hat(ind,:)));
