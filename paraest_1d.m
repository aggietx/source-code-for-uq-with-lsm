%%% parameter estimation 
%%%% the only input is the tracer trajectory
%%%% iterative parameter estimation for impressible flow with reduced modes
gap=1;
disp('paraest_1d')
Nit=6;maxlag=5000;funfit=0;tfit=1;fprintf('lag si %d \n',maxlag);fprintf('total iteration %d \n',Nit);
K_max=Kmax;Kuse=K_max;if K_max~=Kuse;reduced=1;else reduced=0;end;
fprintf('funfit is %d \n',funfit);regularizeval=10^(-10);
estway=2;%%0:filter 1:smoother
constant_cov=0;%%%% use constant covariance matrix.
diagcov=0;if constant_cov==1;diagcov=1;end;if diagcov==0;constant_cov=0;end;
sig_ex1=sig_ex*1.00;%%% use inflation;
dt=dt0;N=N0;x=x(:,gap:gap:end);y=y(:,gap:gap:end);
xexact=x;yexact=y;
%% initial guess
fprintf('Kmax and Kuse are %d %d\n',K_max,Kuse);
if Kuse>K_max;Kuse=K_max;reduced=0;end
dimuhat=(2*K_max+1);

% ix=min(K_max,k_0);iy=min(K_max,k_0);ind=getind(ix,iy,kk);%% sample plot

redind=1:dimuhat;

% u_hat=zeros(1+2*Kmax,N);
% u_hat(1,:)=u;u_hat(2:2:end,:)=psik;
% u_hat(3:2:end,:)=conj(psik);
% [dke,omegake,mue,sigmak2e]=fit_ou_fun(funfit,Kmax,maxlag,dt,u_hat(:,tfit/dt:end),kk1dfullf);
% % load dke;load omegake;load mue;load sigmak2e;
% dkfit=dke;omegakfit=omegake;Ffit=mue;sigmak2=sigmak2e;

dkfit=.1*ones(dimuhat,1);
omegakfit=0.1*ones(dimuhat,1);omegakfit(3:2:end)=-omegakfit(2:2:end);omegakfit(1)=0;
Ffit=0.1*ones(dimuhat,1);
sigmak2=0.1*ones(dimuhat,1);



sigmafit=form_sigmatrix1d(dimuhat,sigmak2,kk1dfullf,K_max);

fprintf('constant_cov and diagcov are %d %d\n',constant_cov,diagcov);
sqrtval=2;%% 2 or 1, 2 represents norm of a+ai,
alldkfit=zeros(dimuhat,1+Nit);
allomegakfit=zeros(dimuhat,1+Nit);
allFfit=zeros(dimuhat,1+Nit);
allsigmak2=zeros(dimuhat,1+Nit);
alldkfit(:,1)=dkfit;allomegakfit(:,1)=omegakfit;
allFfit(:,1)=Ffit;allsigmak2(:,1)=sigmak2;

if constant_cov==1
    disp('constant covariance matrix');
tempsigmafit=sigmak2*sqrt(sqrtval);
fixedcovfit=diag((tempsigmafit.^2)./(dkfit+sqrt(dkfit.^2+L/sig_ex/sig_ex.*tempsigmafit.^2)));
fixedcovfit(1,1)=fixedcovfit(1,1)/sqrt(sqrtval);
fixedcovfit=fixedcovfit(redind,redind);
invRfixed=diag(1./diag(fixedcovfit));invRfixeddiag=1./diag(fixedcovfit);
end
Dim_Y=dimuhat;Dim_Yr=length(redind);
Dim_X=2*L;
kk1dfull=zeros(2,2*K_max);
for i=1:K_max
   kk1dfull(1,2*i-1)=i;
   kk1dfull(1,2*i)=-i;
end
kk1dfullf=[[0;0],kk1dfull];
repk1d=repmat(kk1dfull(1,:),L,1)*1i;
rng(10);
%% parameter estimation
for it=1:Nit
%     close all;
    fprintf('iteration %d\n',it);
Dim_Y=dimuhat;Dim_X=2*L;
gamma_mean0 = zeros(Dim_Yr,1);
% 
if  constant_cov==1
gamma_cov0=fixedcovfit;
else
 gamma_cov0 = eye(Dim_Yr)*0.01;    
end
gamma_mean_trace=zeros(Dim_Y,N);
if estway~=0
    if diagcov
    gamma_cov_trace=zeros(Dim_Yr,N);
    gamma_cov_trace(:,1) = diag(gamma_cov0);        
    else
    gamma_cov_trace=zeros(Dim_Yr,Dim_Yr,N);
    gamma_cov_trace(:,:,1) = gamma_cov0;
    end
end

gamma_mean_trace(redind,1) = gamma_mean0;

A0=zeros(Dim_X,1);
a0=Ffit(redind);a1=diag(-dkfit(redind)+1i*omegakfit(redind));
% a0=zeros(length(redind),1);a1=diag(-dkfit(redind)+0*1i*omegakfit(redind));

b1=sigmafit(redind,redind);
b1b1t=b1 * b1';

invBoB = 1 /sig_ex1/ sig_ex1 * eye(2*L); % inverse of the square of the observational noise
invBoBdiag=1 /sig_ex1/ sig_ex1 * ones(2*L,1);

%% filter
disp('data assimilation (filter)......')
for i = 2:N
    if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
    end
    % observational operator 
    x_loc = [xexact(:,i-1),yexact(:,i-1)];
   
%     G1 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
%     G2 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
% Gy1d = [zeros(L,1),exp(1i * x_loc *  kk1dfull).*repk1d ];
% Gu=[ones(L,1),zeros(L,2*K_max)];
% G1=Gu;G2=(Gu+Gy1d);
[G1,G2]=get1dG1(K_max,x_loc );
% computing the difference between the locations in the Lagrangian
    % tracers; need to consider the cases near the boundaries 
    diff_x1 = xexact(:,i) - xexact(:,i-1); diff_x2 = xexact(:,i) - xexact(:,i-1) + L1; diff_x3 = xexact(:,i) - xexact(:,i-1) - L1;  
    diff_y1 = yexact(:,i) - yexact(:,i-1); diff_y2 = yexact(:,i) - yexact(:,i-1) + L1; diff_y3 = yexact(:,i) - yexact(:,i-1) - L1;  
    diff_xtemp = min(abs(diff_x1), abs(diff_x2)); diff_x_index = min(abs(diff_x3), diff_xtemp);
    diff_ytemp = min(abs(diff_y1), abs(diff_y2)); diff_y_index = min(abs(diff_y3), diff_ytemp);
    diff_x1_index = (diff_x_index == abs(diff_x1)); diff_x2_index = (diff_x_index == abs(diff_x2)); diff_x3_index = (diff_x_index == abs(diff_x3)); 
    diff_y1_index = (diff_y_index == abs(diff_y1)); diff_y2_index = (diff_y_index == abs(diff_y2)); diff_y3_index = (diff_y_index == abs(diff_y3)); 
    diff_x = diff_x1 .* diff_x1_index + diff_x2 .* diff_x2_index + diff_x3 .* diff_x3_index;
    diff_y = diff_y1 .* diff_y1_index + diff_y2 .* diff_y2_index + diff_y3 .* diff_y3_index;
    diff_xy = [diff_x; diff_y];
    
    A1=[G1;G2];A1=A1(:,redind);
    


    % run the data assimilation for posterior mean and posterior covariance
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + ...
        (gamma_cov0 * A1') * (invBoBdiag.* (diff_xy - A0*dt-A1 * gamma_mean0 * dt));
% %     gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + ...
% % %         (gamma_cov0 * A1') * invBoB* (diff_xy - A0*dt-A1 * gamma_mean0 * dt);
    
if constant_cov==1    
gamma_cov=fixedcovfit;%%% use fixed cov matrix
else
%  gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
  gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1b1t - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     

end
    gamma_mean_trace(redind,i) = gamma_mean;
    

    % update
    gamma_mean0 = gamma_mean;    
%     if diagcov
%     gamma_cov0 = diag(diag(gamma_cov));
%     else
    gamma_cov0 = gamma_cov;
%     end
    if diagcov 
    gamma_cov_trace(:,i) = diag(gamma_cov0);
    else
gamma_cov_trace(:,:,i) = gamma_cov0;        
    end
end

if estway==0

[dkfit,omegakfit,Ffit,sigmak2]=fit_ou_fun(funfit,K_max,maxlag,dt,gamma_mean_trace(:,tfit/dt:end),kk1dfullf);
sigmafit=form_sigmatrix1d(dimuhat,sigmak2,kk1dfullf,K_max);
alldkfit(:,1+it)=dkfit;allomegakfit(:,1+it)=omegakfit;
allFfit(:,1+it)=Ffit;allsigmak2(:,1+it)=sigmak2;


    if constant_cov==1
tempsigmafit=sigmak2*sqrt(sqrtval);
fixedcovfit=diag((tempsigmafit.^2)./(dkfit+sqrt(dkfit.^2+L/sig_ex/sig_ex.*tempsigmafit.^2)));
fixedcovfit(1,1)=fixedcovfit(1,1)/sqrt(sqrtval);
fixedcovfit=fixedcovfit(redind,redind);
invRfixed=diag(1./diag(fixedcovfit));invRfixeddiag=1./diag(fixedcovfit);
    end
end
estuhat=gamma_mean_trace(:,1/dt:1/dt:end);
uexact=u(:,1/dt*gap:1/dt*gap:end);
psikexact=psik(:,1/dt*gap:1/dt*gap:end);
rmscc(estuhat(1,:),uexact,1);
rmscc(estuhat(2,:),psikexact(1,:),1);
rmscc(estuhat(4,:),psikexact(2,:),1);
compute_errorlayer1d
errest(:,it)=[mean(rmstime(10:end)),mean(rmsall),mean(ccall)];
 if it==Nit
    break
 end
% return
%% smoother
if estway~=0
   disp('smoother....')
gamma_mean0_smoother = gamma_mean0;
gamma_cov0_smoother = gamma_cov0;
gamma_mean_trace_smoother = zeros(Dim_Y,N);
% gamma_mean_trace_smoother(redind,end)=gamma_mean0;
gamma_mean_trace_smoother(redind,1)=gamma_mean0;
gamma_mean_trace_red=gamma_mean_trace(redind,:);
for i=N:-1:2
   if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',i);
   end
    j=N-i+2;
    % observational operator 
   if constant_cov==1 
% %        invR=invRfixed;
           gamma_mean_smoother=gamma_mean0_smoother+dt*(-a0-a1*gamma_mean0_smoother+...
        (b1b1t)*invRfixeddiag.*(gamma_mean_trace_red(:,i)-gamma_mean0_smoother));
   else
       if diagcov
 R=gamma_cov_trace(:,i);
%  R=regularize(diag(R),regularizeval);
 invR=diag(R)\speye(Dim_Yr);     
% invR=diag(1./R);
       else
    R=gamma_cov_trace(:,:,i);R=regularize(R,regularizeval);invR=R\speye(Dim_Yr);
       end
        gamma_mean_smoother=gamma_mean0_smoother+dt*(-a0-a1*gamma_mean0_smoother+...
        (b1b1t)*invR*(gamma_mean_trace_red(:,i)-gamma_mean0_smoother));
   end
    

% %     gamma_cov_smoother=gamma_cov0_smoother-dt*( (a1+(b1*b1')*invR)*gamma_cov0_smoother+...
% %         gamma_cov0_smoother*(a1'+b1*b1'*invR)-b1*b1' );
        
%     gamma_mean_trace_smoother(redind,i-1) = gamma_mean_smoother;

gamma_mean_trace_smoother(redind,j) = gamma_mean_smoother;
    gamma_mean0_smoother = gamma_mean_smoother;
% %     gamma_cov0_smoother = gamma_cov_smoother;
end
gamma_mean_trace_smoother=gamma_mean_trace_smoother(:,end:-1:1);
    if estway==1


[dkfit,omegakfit,Ffit,sigmak2]=fit_ou_fun(funfit,K_max,maxlag,dt,gamma_mean_trace(:,tfit/dt:end),kk1dfullf);
sigmafit=form_sigmatrix1d(dimuhat,sigmak2,kk1dfullf,K_max);
alldkfit(:,1+it)=dkfit;allomegakfit(:,1+it)=omegakfit;
allFfit(:,1+it)=Ffit;allsigmak2(:,1+it)=sigmak2;

% estuhat=gamma_mean_trace_smoother(:,1/dt:1/dt:end);
% uexact=u(:,1/dt*gap:1/dt*gap:end);
% psikexact=psik(:,1/dt*gap:1/dt*gap:end);
% rmscc(estuhat(1,:),uexact,1);
% rmscc(estuhat(2,:),psikexact(1,:),1);
% rmscc(estuhat(4,:),psikexact(2,:),1);
% compute_errorlayer1d
%  errest(:,it)=[mean(rmstime(10:end)),mean(rmsall),mean(ccall)];

  if constant_cov==1
tempsigmafit=sigmak2*sqrt(sqrtval);
fixedcovfit=diag((tempsigmafit.^2)./(dkfit+sqrt(dkfit.^2+L/sig_ex/sig_ex.*tempsigmafit.^2)));
fixedcovfit(1,1)=fixedcovfit(1,1)/sqrt(sqrtval);
fixedcovfit=fixedcovfit(redind,redind);
invRfixed=diag(1./diag(fixedcovfit));invRfixeddiag=1./diag(fixedcovfit);
 end  
    end
  %% backward sampling
    if estway==2
disp('backward sampling...')
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
     if constant_cov==1 
         invR=invRfixed;
     else
       if diagcov
 R=gamma_cov_trace(:,i);
%  R=regularize(diag(R),regularizeval);
  invR=diag(R)\speye(Dim_Yr);          
       else
    R=gamma_cov_trace(:,:,i);R=regularize(R,regularizeval);invR=R\speye(Dim_Yr);
       end      
     end
    gamma_mean_smoother_sample=gamma_mean0_smoother_sample+dt*(-a0-a1*gamma_mean0_smoother_sample+...
        (b1b1t)*invR*(gamma_mean_trace_red(:,i)-gamma_mean0_smoother_sample))+sqrt(dt)*b1*randn(Dim_Yr,1);
% gamma_mean_trace_smoother_sample(redind,i-1) = gamma_mean_smoother_sample;
gamma_mean_trace_smoother_sample(redind,j) = gamma_mean_smoother_sample;
    gamma_mean0_smoother_sample = gamma_mean_smoother_sample;
end
   
 gamma_mean_trace_smoother_sample=gamma_mean_trace_smoother_sample(:,end:-1:1);    
% ix=k_0;iy=k_0;ind=getind(ix,iy,kk);
% rnorm(real(gamma_mean_trace_smoother_sample(ind,:)),real(u_hat(ind,:)));
[dkfit,omegakfit,Ffit,sigmak2]=fit_ou_fun(funfit,K_max,maxlag,dt,gamma_mean_trace_smoother_sample(:,tfit/dt:end),kk1dfullf);
sigmafit=form_sigmatrix1d(dimuhat,sigmak2,kk1dfullf,K_max);
alldkfit(:,1+it)=dkfit;allomegakfit(:,1+it)=omegakfit;
allFfit(:,1+it)=Ffit;allsigmak2(:,1+it)=sigmak2;

% [dkfit,omegakfit,Ffit,sigmak2]=fit_ou_fun(funfit,1,maxlag,dt,traj(:,1/dt:end),kk1dfullf);

% all_bsmean(:,:,it)=gamma_mean_trace_smoother_sample;
if constant_cov==1
tempsigmafit=sigmak2*sqrt(sqrtval);
fixedcovfit=diag((tempsigmafit.^2)./(dkfit+sqrt(dkfit.^2+L/sig_ex/sig_ex.*tempsigmafit.^2)));
fixedcovfit(1,1)=1;
fixedcovfit=fixedcovfit(redind,redind);
invRfixed=diag(1./diag(fixedcovfit));invRfixeddiag=1./diag(fixedcovfit);
end
%  [rmstimeest,rmsallest,ccallest ]=compute_velerrorGBsw(u_hat(:,gap:gap:end),gamma_mean_trace_smoother(:,gap:gap:end),rk,kk,L1,Dim_Grid,3);
%  errest(:,it)=[mean(rmstimeest(3,:)),mean(rmsallest),mean(ccallest)];

    end%%% end of use backward sampling estimation  
end

% p10=p*10;
% fname = sprintf('rmstime_est1D_it%dp%dL%d.bin',it,p10,L);
% fid=fopen(fname,'w+');fwrite(fid,rmstime,'single');fclose(fid);%%%

end%%%% end of iterative estimation
uhate=zeros(1+2*Kmax,T);
uhate(1,:)=uexact;
uhate(2:2:end,:)=psikexact;
uhate(3:2:end,:)=conj(psikexact);
save estuhat estuhat;save uhate uhate;
fprintf('test case is %d\n',testcase);
fprintf('lag si %d \n',maxlag);fprintf('total iteration %d \n',Nit);
fprintf('Kmax and Kuse are %d %d\n',K_max,Kuse);
fprintf('funfit is %d \n',funfit);
fprintf('L is %d \n',L);
fprintf('constant_cov and diagcov are %d %d\n',constant_cov,diagcov);
eneest1d=sum((real(estuhat)).^2+(imag(estuhat)).^2,2);



% estuhat=gamma_mean_trace_smoother(:,1/dt:1/dt:end);
% uexact=u(:,1/dt:1/dt:end);
% psikexact=psik(:,1/dt:1/dt:end);
% return

% clear x y xexact yexact u psik gamma_cov_trace  gamma_mean_trace 
% clear gamma_mean_trace_red gamma_mean_trace_smoother
% clear fixedcovfit gamma_cov gamma_cov0   tempsigmafit sigmafit
% clear gamma_cov0_smoother invRfixed b1b1t Ga Gb A1 invBoB  invRfixeddiag
% clear dkfit G1 G2 gamma_mean_trace_smoother_sample
% % compute_errorlayer1d
% save data

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
return
[dkfit,omegakfit,Ffit,sigmak2]=fit_ou_fun(funfit,K_max,maxlag,dt,gamma_mean_trace_smoother_sample(:,tfit/dt:end),kk1dfullf);
sigmafit=form_sigmatrix1d(dimuhat,sigmak2,kk1dfullf,K_max);
traj1=zeros(1+2*Kmax,N);traj1(1,:)=u;
traj1(2:2:end,:)=psik;traj1(3:2:end,:)=conj(psik);
[dkfit,omegakfit,Ffit,sigmak2]=fit_ou_fun(funfit,K_max,maxlag,dt,traj1(:,tfit/dt:end),kk1dfullf);
sigmafit=form_sigmatrix1d(dimuhat,sigmak2,kk1dfullf,K_max);


ufree=zeros(dimuhat,1);temp=ufree(:,1);
for i=2:N
    temp1=temp+(-dkfit+1i*omegakfit).*temp*dt+Ffit*dt+sqrt(dt)*sigmafit*randn(dimuhat,1);
    temp=temp1;
    ufree(:,i)=temp;
end
id=2;
[prob,r1]=plotpdf(real(ufree(id,1/dt:end)));a=[prob;r1];save a a
%%%%%%%%%%%%%%%%%%%%
dkfita=dkfit([1,2:2:end]);
omegakfita=omegakfit([1,2:2:end]);
Ffita=Ffit([1,2:2:end]);
sigmak2a=sigmak2([1,2:2:end]);
sigmak=zeros(Kmax+1,1);
sigmak(1)=sigmak2a(1);
sigmak(2:end)=sigmak2a(2:end)+1i*sigmak2a(2:end);

ufree=zeros(Kmax+1,1);temp=ufree(:,1);
for i=2:N
    temp1=temp+(-dkfita+1i*omegakfita).*temp*dt+Ffita*dt+sqrt(dt)*sigmak.*randn(Kmax+1,1);
    temp=temp1;
    ufree(:,i)=temp;
end
id=2;
[prob,r1]=plotpdf(real(ufree(id,1/dt:end)));a1=[prob;r1];save a1 a1