   
%% etkf
bdryh=1;
% % temp_obs=[xexact(:,i-1);yexact(:,i-1)]+0.1*sigma_obs*randn(2*L,1);
temp_obs=[xexact(:,i);yexact(:,i)];
x_pri=[xexact(:,i-1);yexact(:,i-1)]+(A0+A1*gamma_mean0)*dt;
x_pri=mvnrnd(x_pri,sig_ex*speye(2*L)*dt,MM);x_pri=x_pri.';x_pri=real(x_pri);
xensnew=x_pri(1:L,:);yensnew=x_pri(1+L:2*L,:);
mu_pri=gamma_mean0+(a0 + a1 * gamma_mean0) * dt;
R_pri=gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1b1t )*dt;
% % % u_hatpri = mvnrnd(mu_pri,R_pri,MM);u_hatpri=u_hatpri';%%% generate pri ensemble
 u_hatpri1 = mvnrnd(mu_pri(halfind),R_pri(halfind,halfind),MM);u_hatpri1=u_hatpri1.';
u_hatpri=zeros(dimuhat,MM);
u_hatpri(halfind,:)=u_hatpri1;u_hatpri(halfind1,:)=conj(u_hatpri1);
uhatmesh=[invG1;invG2;invG3]*u_hatpri;%%%% back to fourier domain
if max(imag(abs(uhatmesh(:))))>10^(-8)
    disp('ensemble not symmetric')
end
uhatnew=uhatmesh;
   [avexens,xensnew]=compute_ensave_bd(xensnew,bdryh,MM);
  [aveyens,yensnew]=compute_ensave_bd(yensnew,bdryh,MM);
idbdobsx=find(abs(avexens-temp_obs(1:L))>2*pi-2*bdryh);
idbdobsy=find(abs(aveyens-temp_obs(1+L:end))>2*pi-2*bdryh);
temp_obs(idbdobsx)=temp_obs(idbdobsx)-2*pi;
temp_obs(idbdobsy+L)=temp_obs(idbdobsy+L)-2*pi;

   u_prior=[xensnew;yensnew;uhatnew];
    u_mean_prior=sum(u_prior,2)/MM;
    U=u_prior - u_mean_prior * ones(1, MM);
    U=sqrt(1+r)*U;
    V=U(1:2*L,:);%%%% observation residual 
pricov(:,:,i)=cov(U(2*L+1:end,:).');


J = (MM - 1) / (1 + r) * eye(MM) + V' / Ro * V;   
J = (J + J')/2; 
[X, Gamma] = eig(J);   
xx = J \ V' / Ro * ( temp_obs - H * mean(u_prior,2) );
    u_mean_posterior = u_mean_prior + U * xx; % posterior mean
% u_mean_posterior
    T1 = sqrt(MM-1) * X * Gamma^(-1/2) * X';% transform matrix
    U_perturb_posterior = U * T1; % posterior perturbation matrix
    u_posterior = u_mean_posterior * ones(1, MM) + U_perturb_posterior; % posterior ensembles
Upost=u_posterior(1+2*L:end,:);
    covph1=cov(Upost.');
postcov(:,:,i)=covph1;