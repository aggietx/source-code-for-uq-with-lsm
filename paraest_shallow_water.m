%%%% iterative parameter estimation for shallow water with reduced modes
gap1=gap;%%% gap of save data
Nit=8;maxlag1=5000;maxlag2=1000;tfit=1;funfit=1;%%%% actually we use mixed fit (fun(gr) and not fun(gb))
fprintf('lag are %d %d\n',maxlag1,maxlag2);fprintf('Nit is %d\n',Nit);
constant_cov=0;%%%% use constant covariance matrix.
diagcov=0;if constant_cov==1;diagcov=1;end;if diagcov==0;constant_cov=0;end;
estway=1;%%0:filter 1:smoother 2:backward
initcase=1;%%% 1:guess, 2: exact uhat fit, 3: exact parameter smoother fit,4:exact parameter
sig_ex1=sig_ex*1.00;%%% use inflation;
% % [dke,omegake,mue,sigmak2e]=fit_ou_fun_sw(funfit,K_max,maxlag1,maxlag2,dt,u_hat(:,tfit/dt:end),kk);
% alldkfit=[alldkfit,dke];allsigmak2=[allsigmak2,sigmak2e];
% ix=min(K_max,k_0);iy=min(K_max,k_0);ind=getind(ix,iy,kk);%% sample plot
ix=K_max;iy=K_max;ind=getind(ix,iy,kk);%% sample plot
fprintf('initcase is %d\n',initcase);
fprintf('L is %d\n',L);
%% initial guess
if initcase==1
close all;
dkfit=.1*ones(size(u_hat,1),1);dkfit([1,1+dimuhat0,1+2*dimuhat0])=0;
omegakfit0=0.1;omegakfitr=formomega(dimuhat0,omegakfit0,kk0,Kmax);
omegakfit=[omegakfitr;omegakfit0*ones(dimuhat0,1);-omegakfit0*ones(dimuhat0,1)];
omegakfit([1,1+dimuhat0,1+2*dimuhat0])=0;
Ffit=0.0*ones(dimuhat,1);Ffit([1,1+dimuhat0,1+2*dimuhat0])=0;
sigmak2=0.1*ones(dimuhat0,1);sigmak2(1)=0;
sigmafit=zeros(dimuhat,dimuhat);
sigmafit(1:dimuhat0,1:dimuhat0)=form_sigmatrix1(dimuhat0,sigmak2,kk0,Kmax);%%%% already divided by sqrt(2)
sigmafit(1+dimuhat0:end,1+dimuhat0:end)=form_gravity_sigma(kk0,sqrt(2)*sigmak2);
sigma2eb=sigmak2;sigma2eg=sigmak2;

elseif initcase==2
%%%% use parameter with exact uhat fit
[dke,omegake,mue,sigmak2e]=fit_ou_fun_sw(funfit,K_max,maxlag1,maxlag2,dt,u_hat(:,tfit/dt:end),kk);
dkfit=dke;omegakfit=omegake;Ffit=mue;
sigma2eb=sigmak2e(1:dimuhat0);
sigma2eg=sigmak2e(1+dimuhat0:2*dimuhat0);
sigmafit=zeros(dimuhat,dimuhat);
sigmafit(1:dimuhat0,1:dimuhat0)=form_sigmatrix1(dimuhat0,sigma2eb,kk0,Kmax);%%%% already divided by sqrt(2)
sigmafit(1+dimuhat0:end,1+dimuhat0:end)=form_gravity_sigma(kk0,sqrt(2)*sigma2eg);
Nit=1;
elseif initcase==3
    %%%% use exact parameter smoother
[dke,omegake,mue,sigmak2e]=fit_ou_fun_sw(funfit,K_max,maxlag1,maxlag2,dt,exact_gamma_mean_trace(:,tfit/dt:end),kk);
dkfit=dke;omegakfit=omegake;Ffit=mue;
sigma2eb=sigmak2e(1:dimuhat0);
sigma2eg=sigmak2e(1+dimuhat0:2*dimuhat0);
sigmafit=zeros(dimuhat,dimuhat);
sigmafit(1:dimuhat0,1:dimuhat0)=form_sigmatrix1(dimuhat0,sigma2eb,kk0,Kmax);%%%% already divided by sqrt(2)
sigmafit(1+dimuhat0:end,1+dimuhat0:end)=form_gravity_sigma(kk0,sqrt(2)*sigma2eg);
Nit=1;
else
%%% use exact parameter
dkfit=-dkexact;omegakfit=omegaexact;Ffit=Fuexact;
sigmafit=sigmauhat;sigma2eb=sigmak2exactb/sqrt(2);sigma2eg=sigmak2exactg/sqrt(2);
Nit=1;
end
fprintf('Nit is %d\n',Nit);


fprintf('constanat cov, diag cov: %d %d\n',constant_cov,diagcov);
sqrtval=2;%% 2 or 1, 2 represents norm of a+ai,
if constant_cov==1
    disp('constant covariance matrix');
% tempsigmafit=[sigmak2;sigmak2;sigmak2]*sqrt(sqrtval);
tempsigmafit=[sigma2eb;sigma2eg;sigma2eg]*sqrt(sqrtval);
fixedcovfit=diag((tempsigmafit.^2)./(dkfit+sqrt(dkfit.^2+L/sig_ex/sig_ex.*tempsigmafit.^2)));
fixedcovfit(1,1)=1;fixedcovfit(1+dimuhat0,1+dimuhat0)=1;fixedcovfit(1+2*dimuhat0,1+2*dimuhat0)=1;
fixedcovfit=fixedcovfit(redind,redind);
invRfixed=diag(1./diag(fixedcovfit));invRfixeddiag=1./diag(fixedcovfit);
end
% return
%% record data
alldkfit=zeros(dimuhat,2+Nit);
allomegakfit=zeros(dimuhat,2+Nit);
allFfit=zeros(dimuhat,2+Nit);
allsigmak2=zeros(dimuhat,2+Nit);
alldkfit(:,1)=dkfit;allomegakfit(:,1)=omegakfit;allFfit(:,1)=Ffit;
allsigmak2(:,1)=[sigma2eb;sigma2eg;sigma2eg];
alldkfit(:,end)=dkexact;allomegakfit(:,end)=omegaexact;allFfit(:,end)=Fuexact;
allsigmak2(:,end)=[sigmak2exactb;sigmak2exactg;sigmak2exactg]/sqrt(2);

all_filmean=zeros(dimuhat,N,1+Nit);all_filmean(:,:,end)=gamma_mean_trace;
if exist('gamma_mean_trace_smoother','var') && estway==1 
all_smmean=zeros(dimuhat,N,1+Nit);all_smmean(:,:,end)=gamma_mean_trace_smoother;
end
if  exist('gamma_mean_trace_smoother_sample','var') && estway==2
 all_bsmean=zeros(dimuhat,N,1+Nit);all_bsmean(:,:,end)=gamma_mean_trace_smoother_sample;  
end
% % utest=zeros(dimuhat,100);%%%% test conjuate
% % for i=3:100
% %        utest(:,i)=utest(:,i-1)+(dkfit+1i*omegakfit).*utest(:,i-1)*dt+...
% %         Ffit*dt+sqrt(dt)*sigmafit*randn(dimuhat,1); 
% %     x_loc = [xexact(:,i-1),yexact(:,i-1)];
% %     G1 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
% %     G2 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
% %     u = (G1*  utest(:,i-1)); 
% %     v = (G2*  utest(:,i-1)); 
% %     if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
% %       disp('complex velocity'  )
% %       max(abs(imag(u)))
% %     end  
% % end
% % 
% % return
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
a0=Ffit(redind);
a1=diag(-dkfit(redind)+1i*omegakfit(redind));
b1=sigmafit(redind,redind);

% % A0=zeros(Dim_X,1);
% % a0=Fuexact(redind);
% % a1=diag(dkexact(redind)+1i*omegaexact(redind));
% % b1=sigmauhat(redind,redind);


b1b1t=b1 * b1';
x=xexact;y=yexact;
invBoB = 1 /sig_ex1/ sig_ex1 * eye(2*L); % inverse of the square of the observational noise
invBoBdiag=1 /sig_ex1/ sig_ex1 * ones(2*L,1);
time=1;
%% filter
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
% ix=k_0;iy=k_0;ind=getind(ix,iy,kk);ind0=ind(1);
% rnorm(real(gamma_mean_trace(ind0,:)),real(u_hat(ind0,:)));
rmscc(real(gamma_mean_trace(ind(1),:)),real(u_hat(ind(1),:)),1);
rmscc(real(gamma_mean_trace(ind(2),:)),real(u_hat(ind(2),:)),1);


[dkfit,omegakfit,Ffit,sigmak2]=fit_ou_fun_sw(funfit,K_max,maxlag1,maxlag2,dt,gamma_mean_trace(:,tfit/dt:end),kk);
sigmafit=zeros(dimuhat,dimuhat);
sigmafit(1:dimuhat0,1:dimuhat0)=form_sigmatrix1(dimuhat0,sigmak2(1:dimuhat0),kk0,Kmax);
sigmafit(1+dimuhat0:end,1+dimuhat0:end)=form_gravity_sigma(kk0,sqrt(2)*sigmak2(1+dimuhat0:2*dimuhat0));
% norm( sigmak2(negindkk+2*dimuhat0)-sigmak2(dimuhat0+1:2*dimuhat0))

alldkfit(:,it+1)=dkfit;allomegakfit(:,it+1)=omegakfit;allFfit(:,it+1)=Ffit;
allsigmak2(:,it+1)=sigmak2;
all_filmean(:,:,it)=gamma_mean_trace;
    if constant_cov==1
tempsigmafit=sigmak2*sqrt(sqrtval);
fixedcovfit=diag((tempsigmafit.^2)./(dkfit+sqrt(dkfit.^2+L/sig_ex/sig_ex.*tempsigmafit.^2)));
fixedcovfit(1,1)=1;fixedcovfit(1+dimuhat0,1+dimuhat0)=1;fixedcovfit(1+2*dimuhat0,1+2*dimuhat0)=1;
fixedcovfit=fixedcovfit(redind,redind);
invRfixed=diag(1./diag(fixedcovfit));invRfixeddiag=1./diag(fixedcovfit);
    end
end


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
% ix=k_0;iy=k_0;ind=getind(ix,iy,kk);
% rnorm(real(gamma_mean_trace_smoother(ind,:)),real(u_hat(ind,:)));
rmscc(real(gamma_mean_trace_smoother(ind(1),:)),real(u_hat(ind(1),:)),1);
rmscc(real(gamma_mean_trace_smoother(ind(2),:)),real(u_hat(ind(2),:)),1);

[dkfit,omegakfit,Ffit,sigmak2]=fit_ou_fun_sw(funfit,K_max,maxlag1,maxlag2,dt,gamma_mean_trace_smoother(:,tfit/dt:end),kk);
sigmafit=zeros(dimuhat,dimuhat);
sigmafit(1:dimuhat0,1:dimuhat0)=form_sigmatrix1(dimuhat0,sigmak2(1:dimuhat0),kk0,Kmax);
sigmafit(1+dimuhat0:end,1+dimuhat0:end)=form_gravity_sigma(kk0,sqrt(2)*sigmak2(1+dimuhat0:2*dimuhat0));



alldkfit(:,it+1)=dkfit;allomegakfit(:,it+1)=omegakfit;allFfit(:,it+1)=Ffit;
allsigmak2(:,it+1)=sigmak2;   
all_smmean(:,:,it)=gamma_mean_trace_smoother;
  if constant_cov==1
tempsigmafit=sigmak2*sqrt(sqrtval);
fixedcovfit=diag((tempsigmafit.^2)./(dkfit+sqrt(dkfit.^2+L/sig_ex/sig_ex.*tempsigmafit.^2)));
fixedcovfit(1,1)=1;fixedcovfit(1+dimuhat0,1+dimuhat0)=1;fixedcovfit(1+2*dimuhat0,1+2*dimuhat0)=1;
fixedcovfit=fixedcovfit(redind,redind);
invRfixed=diag(1./diag(fixedcovfit));invRfixeddiag=1./diag(fixedcovfit);
  end  
 [rmstimesm,rmsallsm,ccallsm ]=compute_velerrorGBsw(u_hat(:,gap:gap:end),gamma_mean_trace_smoother(:,gap:gap:end),rk,kk,L1,Dim_Grid,3);
 errest(:,it)=[mean(rmstimesm(3,:)),mean(rmsallsm),mean(ccallsm)];
 
    end

end
end%%%% end of iterative estimation
disp('...........................')
disp('exact parameter smoother sw error');
rmscc(real(exact_gamma_mean_trace_smoother(ind(1),:)),real(u_hat(ind(1),:)),1);
rmscc(real(exact_gamma_mean_trace_smoother(ind(2),:)),real(u_hat(ind(2),:)),1);% rnorm(real(gamma_mean_trace(ind0,:)),real(u_hat(ind0,:)));

 [rmstimeexactsm,rmsexactsm,ccexactsm  ]=compute_velerrorGBsw(u_hat(:,gap:gap:end),exact_gamma_mean_trace_smoother(:,gap:gap:end),rk,kk,L1,Dim_Grid,3);
 if diagcov~=1
 [sigsm,dissm]=computeRentropy(u_hat(:,gap:gap:end),exact_gamma_mean_trace_smoother(:,gap:gap:end),exact_gamma_cov_trace(:,:,gap:gap:end),Rexact,diagcov);
 else
 [sigsm,dissm]=computeRentropy(u_hat(:,gap:gap:end),exact_gamma_mean_trace_smoother(:,gap:gap:end),exact_gamma_cov_trace(:,gap:gap:end),Rexact,diagcov);    
 end
 errsm=[mean(rmstimeexactsm(3:end));mean(rmsexactsm);mean(ccexactsm)];

disp('...........................')
disp('estimated parameter smoother sw error');
fprintf('constanat cov, diag cov: %d %d\n',constant_cov,diagcov);
 [rmstimeest,rmsallest,ccallest ]=compute_velerrorGBsw(u_hat(:,gap:gap:end),gamma_mean_trace_smoother(:,gap:gap:end),rk,kk,L1,Dim_Grid,3);
 if diagcov~=1
 [sigest,disest]=computeRentropy(u_hat(:,gap:gap:end),gamma_mean_trace_smoother(:,gap:gap:end),gamma_cov_trace(:,:,gap:gap:end),Rexact,diagcov);
 else
 [sigest,disest]=computeRentropy(u_hat(:,gap:gap:end),gamma_mean_trace_smoother(:,gap:gap:end),gamma_cov_trace(:,gap:gap:end),Rexact,diagcov);    
 end
% if estway==1
% smest=real(gamma_mean_trace_smoother(ind(2),:));
% ex=real(u_hat(ind(2),:));
% save smest smest;save ex ex
% end
[dke,omegake,mue,sigmak2e]=fit_ou_fun_sw(funfit,K_max,maxlag1,maxlag2,dt,u_hat(:,tfit/dt:end),kk);
alldkfit=[alldkfit,dke];allomegakfit=[allomegakfit,omegake];
allFfit=[allFfit,mue];allsigmak2=[allsigmak2,sigmak2e];
save alldkfit alldkfit;save allomegakfit allomegakfit; save allsigmak2 allsigmak2
fprintf('lag are %d %d\n',maxlag1,maxlag2);
% disp('---------------------------------------------------------------------------------');
eL2dk=zeros(Nit+1,1);eL2omega=zeros(Nit+1,1);eL2sigma=zeros(Nit+1,1);
eL2F=zeros(Nit+1,1);
for ii=1:Nit+1
   ex=-alldkfit(redind,end-1); 
   app=alldkfit(redind,ii);
   eL2dk(ii)=norm(app-ex)/norm(ex);
 
    ex=allomegakfit(redind,end-1); 
   app=allomegakfit(redind,ii);
   eL2omega(ii)=norm(app-ex)/norm(ex);

   app=allFfit(redind,ii);
eL2F(ii)=norm(app);

   ex=allsigmak2(redind,end-1); 
   app=allsigmak2(redind,ii);
   eL2sigma(ii)=norm(app-ex)/norm(ex); 
   
end
fprintf('L2 error of dk is %2.3f\n',eL2dk);
fprintf('L2 error of omegak is %2.3f\n',eL2omega);
fprintf('L2 error of sigma is %2.3f\n',eL2sigma);

eL2dke=zeros(1,1);eL2omegae=zeros(1,1);eL2sigmae=zeros(1,1);
for ii=1:1
   ex=-alldkfit(redind,end-1); 
   app=alldkfit(redind,end);
   eL2dke(ii)=norm(app-ex)/norm(ex);
 
    ex=allomegakfit(redind,end-1); 
   app=allomegakfit(redind,end);
   eL2omegae(ii)=norm(app-ex)/norm(ex);
   
   ex=allsigmak2(redind,end-1); 
   app=allsigmak2(redind,end);
   eL2sigmae(ii)=norm(app-ex)/norm(ex); 
   
end
fprintf('exact fit L2 error of dk is %2.3f\n',eL2dke);
fprintf('exact fit L2 error of omegak is %2.3f\n',eL2omegae);
fprintf('exact fit L2 error of sigma is %2.3f\n',eL2sigmae);
fprintf('L is %d\n',L);
% gap1=10;
savedtraj=[u_hat(:,gap1:gap1:end);
    exact_gamma_mean_trace_smoother(:,1:gap1:end);
    gamma_mean_trace_smoother(:,1:gap1:end)];
compute_freerunsw
disp('---------------------------------------------------------------------------------');

return
%% plot estimation comparison
figure();
plot(-dkexact(redind),'*-r','linewidth',2);hold on;
plot(alldkfit(redind,1),'*-b','linewidth',2);hold on;
plot(alldkfit(redind,end),'*-g','linewidth',2);
lgnd=legend('Exact','Initial','Estimated');
set(lgnd,'FontSize',18);
setgca(18)
title('d_k','FontSize',18);

% figure();
% plot(real(sigmak2exact(redind))/sqrt(2),'*-r','linewidth',2);hold on;
% plot(real(allsigmak2((redind),1)),'*-b','linewidth',2);
% plot(real(allsigmak2((redind),end)),'*-g','linewidth',2);
% lgnd=legend('Exact','Initial','Estimated');
% set(lgnd,'FontSize',18);
% setgca(18)
% title('\sigma_k','FontSize',18);
return
figure();
plot(omegaexact,'*-r','linewidth',2);hold on;
plot(omegakfit,'*-g','linewidth',2);
lgnd=legend('Exact','Estimated');
set(lgnd,'FontSize',18);
setgca(18)
title('\omega_k','FontSize',18);


figure();
plot(real(Fuexact),'*-r','linewidth',2);hold on;
plot(real(Ffit),'*-g','linewidth',2);
lgnd=legend('Exact','Estimated');
set(lgnd,'FontSize',18);
setgca(18)
title('F_k','FontSize',18);

return











% 
%  %% backward sampling
%     if estway==2
% disp('backward sampling...')
% gamma_mean0_smoother_sample = gamma_mean0;
% gamma_mean_trace_smoother_sample = zeros(Dim_Y,N);
% % gamma_mean_trace_smoother_sample(redind,end)=gamma_mean0;
% gamma_mean_trace_smoother_sample(redind,1)=gamma_mean0;
% for i=N:-1:2
%    if mod(i,floor(N/5)) == 0
% %         disp(i*dt)
%  fprintf('rest step is %d\n',i);
%    end  
%    j=N-i+2;
%      if constant_cov==1 
%          invR=invRfixed;
%      else
%        if diagcov
%  R=gamma_cov_trace(:,i);
% %  R=regularize(diag(R),regularizeval);
%  invR=R\speye(Dim_Yr);           
%        else
%     R=gamma_cov_trace(:,:,i);R=regularize(R,regularizeval);invR=R\speye(Dim_Yr);
%        end      
%      end
%     gamma_mean_smoother_sample=gamma_mean0_smoother_sample+dt*(-a0-a1*gamma_mean0_smoother_sample+...
%         (b1b1t)*invR*(gamma_mean_trace_red(:,i)-gamma_mean0_smoother_sample))+sqrt(dt)*b1*randn(Dim_Yr,1);
% % gamma_mean_trace_smoother_sample(redind,i-1) = gamma_mean_smoother_sample;
% gamma_mean_trace_smoother_sample(redind,j) = gamma_mean_smoother_sample;
%     gamma_mean0_smoother_sample = gamma_mean_smoother_sample;
% end
%    
%  gamma_mean_trace_smoother_sample=gamma_mean_trace_smoother_sample(:,end:-1:1);    
% ix=k_0;iy=k_0;ind=getind(ix,iy,kk);
% rnorm(real(gamma_mean_trace_smoother_sample(ind,:)),real(u_hat(ind,:)));
% [dkfit,omegakfit,Ffit,sigmak2]=fit_ou_fun_sw(funfit,K_max,maxlag1,maxlag2,dt,gamma_mean_trace_smoother_sample(:,tfit/dt:end),kk);
% sigmafit=form_sigmatrix1(size(u_hat,1),sigmak2,kk,Kmax);
% alldkfit(:,it+1)=dkfit;allomegakfit(:,it+1)=omegakfit;allFfit(:,it+1)=Ffit;
% allsigmak2(:,it+1)=sigmak2;  
% all_bsmean(:,:,it)=gamma_mean_trace_smoother_sample;
% if constant_cov==1
% tempsigmafit=sigmak2*sqrt(sqrtval);
% fixedcovfit=diag((tempsigmafit.^2)./(dkfit+sqrt(dkfit.^2+L/sig_ex/sig_ex.*tempsigmafit.^2)));
% fixedcovfit(1,1)=1;
% fixedcovfit=fixedcovfit(redind,redind);
% invRfixed=diag(1./diag(fixedcovfit));invRfixeddiag=1./diag(fixedcovfit);
% end
%     end%%% end of use backward sampling estimation   
