MM=200;sigma_obs=0.1;sig_en=0.1;r=0;
fprintf('MM is %d\n',MM);
fprintf('L is %d\n',L);
Dt=1*dt;gap=Dt/dt;
fprintf('Dt and dt are %2.2e %2.2e\n',Dt,dt);
bdryh=1;omega1=H1/sqrt(2);omega3=H2/sqrt(2);
    du=0.0125;dv1=0.0125;dv2=0.0125;dv3=0.0125;dv4=0.0125;

uens=zeros(1,MM);
v1ens=zeros(1,MM);
v2ens=zeros(1,MM);
v3ens=zeros(1,MM);
v4ens=zeros(1,MM);
xobs=x+sigma_obs/2/pi*L1*randn(L,N);
yobs=y+sigma_obs/2/pi*L1*randn(L,N);
xens=repmat(x(:,1),1,MM)+sig_en/2/pi*L1*randn(L,MM);
yens=repmat(y(:,1),1,MM)+sig_en/2/pi*L1*randn(L,MM);
Ro=(sigma_obs*L1/2/pi).^2*speye(2*L);
H = eye(2*L,2*L+5);         % observation matrix
uhatmean=zeros(5,N/gap);
G1=zeros(5,5);%%% (u,v1,v2,v3,v4) to (u,pis1 pis-1,psi2,psi-2);
G1(1,1)=1;
G1(2,2)=1/2/sqrt(2)*(-1-1i);G1(2,3)=1/2/sqrt(2)*(1-1i);
G1(3,2)=1/2/sqrt(2)*(-1+1i);G1(3,3)=1/2/sqrt(2)*(1+1i);
G1(4,4)=1/2/sqrt(2)*(-1-1i);G1(4,5)=1/2/sqrt(2)*(1-1i);
G1(5,4)=1/2/sqrt(2)*(-1+1i);G1(5,5)=1/2/sqrt(2)*(1+1i);
    kk1dfull=[1,-1,2,-2;0,0,0,0];
    repk1d=repmat([1,-1,2,-2],L,1)*1i;
for ii=2:N
        if mod(ii,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-ii);
        end   
uensnew=uens+(omega1*v1ens+2*omega3*v3ens-du*uens)*dt+...
    sigmau*randn(1,MM)*sqrt(dt);

v1ensnew=v1ens+(-beta*v2ens+v2ens.*uens-...
    2*omega1*uens-dv1*v1ens)*dt+sigmau*randn(1,MM)*sqrt(dt);

v2ensnew=v2ens+(beta*v1ens+v1ens.*uens-...
    dv2*v2ens)*dt+sigmau*randn(1,MM)*sqrt(dt);  

v3ensnew=v3ens+(-beta/2*v4ens+2*v4ens.*uens-...
     omega3*uens-dv3*v3ens)*dt+sigmau*randn(1,MM)*sqrt(dt); 
 
v4ensnew=v4ens+(beta/2*v3ens-2*v3ens.*uens-...
    dv4*v4ens)*dt+sigmau*randn(1,MM)*sqrt(dt); 
ufull=zeros(L,MM);vfull=zeros(L,MM);
for j=1:MM
  x_loc = [xens(:,j),yens(:,j)];%+sig_ex*sqrt(dt)*randn(L,2);  
Gy1d = [zeros(L,1),exp(1i * x_loc *  kk1dfull).*repk1d ];
Gu=[ones(L,1),zeros(L,4)];
Gux=Gu*G1;Guy=(Gu+Gy1d)*G1;
gx1=Gux*[uens(j),v1ens(j),v2ens(j),v3ens(j),v4ens(j)]';
gx2=Guy*[uens(j),v1ens(j),v2ens(j),v3ens(j),v4ens(j)]';
 ufull (:,j)=gx1;
 vfull (:,j)=gx2;
end

% whos uens
% pause()
if max(abs(imag(ufull(:))))>10^(-6) || max(abs(imag(vfull(:))))>10^(-6)
    disp('error, complex vel 1');
else
    ufull=real(ufull);vfull=real(vfull);
end

    xensnew = xens + ufull * dt + sqrt(dt) * sig_ex * randn(L,MM); %%% predicted new x
    yensnew = yens + vfull * dt + sqrt(dt) * sig_ex * randn(L,MM); %%% predicted new y
xensnew=real(xensnew);yensnew=real(yensnew);

xensnew = mod(xensnew,L1);  yensnew = mod(yensnew,L1);

if mod(ii,gap)==0
   [avexens,xensnew]=compute_ensave_bd2(xensnew,bdryh*L1/(2*pi),MM,L1);
  [aveyens,yensnew]=compute_ensave_bd2(yensnew,bdryh*L1/(2*pi),MM,L1);
 temp_obs=[xobs(:,ii);yobs(:,ii)]; 
idbdobsx=find(abs(avexens-temp_obs(1:L))>L1-2*bdryh*L1/(2*pi));
idbdobsy=find(abs(aveyens-temp_obs(1+L:end))>L1-2*bdryh*L1/(2*pi));
temp_obs(idbdobsx)=temp_obs(idbdobsx)-L1;
temp_obs(idbdobsy+L)=temp_obs(idbdobsy+L)-L1;
uhatnew=[uensnew;v1ensnew;v2ensnew;v3ensnew;v4ensnew];
    u_prior=[xensnew;yensnew;uhatnew];
    u_mean_prior=sum(u_prior,2)/MM;
    U=u_prior - u_mean_prior * ones(1, MM);
    U=sqrt(1+r)*U;
    V=H*U;%%%% observation residual

J = (MM - 1) / (1 + r) * eye(MM) + V' / Ro * V;   
J = (J + J')/2; 
[X, Gamma] = eig(J);   
x1 = J \ V' / Ro * ( temp_obs - H * mean(u_prior,2) );
    u_mean_posterior = u_mean_prior + U * x1; % posterior mean
% u_mean_posterior
    T1 = sqrt(MM-1) * X * Gamma^(-1/2) * X';% transform matrix
    U_perturb_posterior = U * T1; % posterior perturbation matrix
    u_posterior = u_mean_posterior * ones(1, MM) + U_perturb_posterior; % posterior ensembles
   v4ens=u_posterior(5+2*L,:); 
   v3ens= u_posterior(4+2*L,:); 
   v2ens= u_posterior(3+2*L,:);
   v1ens= u_posterior(2+2*L,:);
   uens=u_posterior(1+2*L,:);
   xens = u_posterior(1:L,:);
   yens= u_posterior(1+L:2*L,:);  
   uhatmean(:,ii/gap)=mean(u_posterior(2*L+1:end,:),2);
else
   v4ens= v4ensnew; v3ens= v3ensnew; v2ens= v2ensnew; v1ens= v1ensnew;
   uens=uensnew;
   xens = xensnew;
   yens= yensnew;
end
end
psi1mean=1/2/sqrt(2)*( (uhatmean(3,:)-uhatmean(2,:))-1i*(uhatmean(3,:)+uhatmean(2,:)) );
psi2mean=1/2/sqrt(2)*( (uhatmean(5,:)-uhatmean(4,:))-1i*(uhatmean(5,:)+uhatmean(4,:)) );

rmscc(uhatmean(1,:),u(gap:gap:end),1);
rmscc(psi1mean,psik(1,gap:gap:end),1);
rmscc(psi2mean,psik(2,gap:gap:end),1);
ue=u(gap:gap:end);uetkf=uhatmean(1,:);
save ue ue
save uetkf uetkf
% get ue.mat uetkf.mat
% load ue;load uetkf;
% plot2line(ue,uetkf)