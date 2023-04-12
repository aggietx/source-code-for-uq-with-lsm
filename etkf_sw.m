%% etkf initilization
rng(10);MM=200;gamma=20;sigma_obs=0.1;nfloe=L;sigma_init=0.1;
gap=0.1/dt;fprintf('ensemble size is %d\n',MM);
fprintf('gap is %d\n',gap);
disp('etkf without transform, sw')
xobs=xexact+sigma_obs*randn(nfloe,N);%%%% observation
yobs=yexact+sigma_obs*randn(nfloe,N);%%%% observation
uhat=u_hat(:,1)*ones(1,MM)+1*sigma_init*sigmauhat*(randn(dimuhat,MM));%%% initial uhat
xens=xexact(:,1)*ones(1,MM)+sigma_init*(randn(nfloe,MM));%%%% ensemble
yens=yexact(:,1)*ones(1,MM)+sigma_init*(randn(nfloe,MM));%%%% ensemble
repL_u=repmat(dkexact+1i*omegaexact,1,MM);
repFu=repmat(Fuexact,1,MM);
xy=kk0'/(2*K_max+1)*2*pi;
rketa=[1./sqrt(kk0(1,:).^2 + kk0(2,:).^2+1),...
    1./sqrt(kk0(1,:).^2 + kk0(2,:).^2)/sqrt(2)./sqrt(kk0(1,:).^2 + kk0(2,:).^2 + 1).*(kk0(1,:).^2 + kk0(2,:).^2),...
    1./sqrt(kk0(1,:).^2 + kk0(2,:).^2)/sqrt(2)./sqrt(kk0(1,:).^2 + kk0(2,:).^2 + 1).*(kk0(1,:).^2 + kk0(2,:).^2)];
rketa(1+dimuhat0)=0;rketa(1+dimuhat0*2)=0;
    invG1 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rk(1,:))); % Fourier bases for u
    invG2 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rk(2,:))); % Fourier bases for v
    invG3 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rketa)); % Fourier bases for v
r=0;bdryh=1;

    
G=exp(-1i  * xy*kk0)/(2*K_max+1)^2;%%%%% shallow water matrix from physical domain to F domain
    Gm=[G,zeros(dimuhat0,dimuhat0),zeros(dimuhat0,dimuhat0);
        zeros(dimuhat0,dimuhat0),G,zeros(dimuhat0,dimuhat0);
        zeros(dimuhat0,dimuhat0),zeros(dimuhat0,dimuhat0),G];
gamma_mean_trace_local=zeros(dimuhat,N);
Ro=(sigma_obs).^2*speye(2*L);
H=speye(2*L,2*L+dimuhat);
covph=zeros(dimuhat,dimuhat,N/gap);
covphpri=zeros(dimuhat,dimuhat,N/gap);

for i=2:N
%     i
    if mod(i,floor(N/10)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
    end
    %% predict
    uhatnew=uhat+repL_u.*uhat*dt+repFu*dt+1*sqrt(dt)*sigmauhat*randn(dimuhat,MM);
    ufull=zeros(L,MM);vfull=zeros(L,MM);   
    for ii=1:MM
       x_loc=[xens(:,ii),yens(:,ii)]; 
    G1 = (exp(1i * x_loc * kk) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
    G2 = (exp(1i * x_loc * kk) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
    ufull(:,ii) = (G1*  uhat(:,ii)); 
    vfull(:,ii) = (G2*  uhat(:,ii));  
    end
if max(abs(imag(ufull(:))))>10^(-6) || max(abs(imag(vfull(:))))>10^(-6)
    disp('error, complex vel 1');
else
    ufull=real(ufull);vfull=real(vfull);
end
    xensnew = xens + ufull * dt ; % floe equation in x
    yensnew = yens + vfull* dt ; % floe equation in y
    xensnew=mod(xensnew,2*pi);
    yensnew=mod(yensnew,2*pi);
if mod(i,gap)==0
    %% analysis prepare
   [avexens,xensnew]=compute_ensave_bd(xensnew,bdryh,MM);
  [aveyens,yensnew]=compute_ensave_bd(yensnew,bdryh,MM);
 temp_obs=[xobs(:,i);yobs(:,i)]; 
idbdobsx=find(abs(avexens-temp_obs(1:L))>2*pi-2*bdryh);
idbdobsy=find(abs(aveyens-temp_obs(1+L:end))>2*pi-2*bdryh);
temp_obs(idbdobsx)=temp_obs(idbdobsx)-2*pi;
temp_obs(idbdobsy+L)=temp_obs(idbdobsy+L)-2*pi;

   u_prior=[xensnew;yensnew;uhatnew];
    u_mean_prior=sum(u_prior,2)/MM;
    U=u_prior - u_mean_prior * ones(1, MM);
    U=sqrt(1+r)*U;
    V=U(1:2*L,:);%%%% observation residual 
% covphpri(:,:,i/gap)=U(2*L+1:end,:)*U(2*L+1:end,:)';
covphpri(:,:,i/gap)=cov(U(2*L+1:end,:).');

      Ux=U(1:L,:);Uy=U(1+L:2*L,:);
    Uhat=U(1+2*L:end,:);
J = (MM - 1) / (1 + r) * eye(MM) + V' / Ro * V;   
J = (J + J')/2; 
[X, Gamma] = eig(J);   
x = J \ V' / Ro * ( temp_obs - H * mean(u_prior,2) );
    u_mean_posterior = u_mean_prior + U * x; % posterior mean
% u_mean_posterior
    T1 = sqrt(MM-1) * X * Gamma^(-1/2) * X';% transform matrix
    U_perturb_posterior = U * T1; % posterior perturbation matrix
    u_posterior = u_mean_posterior * ones(1, MM) + U_perturb_posterior; % posterior ensembles

     uhat=u_posterior(1+2*L:end,:);
     xens=u_posterior(1:L,:);
     yens=u_posterior(1+L:2*L,:);
    xens=mod(xens,2*pi);
    yens=mod(yens,2*pi);
    gamma_mean_trace_local(:,i)=u_mean_posterior(1+2*L:end);  

% covph(:,:,i/gap)= u_posterior(2*L+1:end,:)* u_posterior(2*L+1:end,:)';
 covph(:,:,i/gap)= cov( u_posterior(2*L+1:end,:).');
else%%% no update
     uhat=uhatnew;
     xens= xensnew;
     yens= yensnew;
    xens=mod(xens,2*pi);
    yens=mod(yens,2*pi);
end
end
% fprintf('gamma is %2.2f\n',gamma);
% % figure();plot(dt:dt:T,real(u_hat(ind(1),:)),'r','linewidth',1.5);hold on;
% % plot(dt:dt:T,real(gamma_mean_trace_local(ind(1),:)),'b','linewidth',1.5);
% % setgca(16);title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
% % xlabel('t','fontsize',16)
rmscc(real(gamma_mean_trace_local(ind(1),gap:gap:end)),real(u_hat(ind(1),gap:gap:end)),1);
rmscc(real(gamma_mean_trace_local(ind(2),gap:gap:end)),real(u_hat(ind(2),gap:gap:end)),1);

% rmscc(real(gamma_mean_trace_local(ind(1),gap:gap:gap*100)),real(u_hat(ind(1),gap:gap:gap*100)),1);
% rmscc(real(gamma_mean_trace_local(ind(2),gap:gap:gap*100)),real(u_hat(ind(2),gap:gap:gap*100)),1);