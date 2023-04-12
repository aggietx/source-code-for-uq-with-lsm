%% letkf initilization should use etkf_gb to generate data first
rng(10)
MM=200;gamma=2.0;sigma_obs=0.1;nfloe=L;sigma_init=0.1;dimuhat0=dimuhat;kk0=kk;
gap=0.1/dt;fprintf('gamma is %2.2f\n',gamma);fprintf('ensemble size is %d\n',MM);
disp('letkf with rk transform,gb')
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
    invG3 = (exp(1i * xy * kk) .* (ones(length(xy),1) * rketa(1:dimuhat0))); % Fourier bases for v
r=0;bdryh=1;

    
G=exp(-1i  * xy*kk0)/(2*K_max+1)^2;%%%%% shallow water matrix from physical domain to F domain
    Gm=[G,zeros(dimuhat0,dimuhat0),zeros(dimuhat0,dimuhat0);
        zeros(dimuhat0,dimuhat0),G,zeros(dimuhat0,dimuhat0);
        zeros(dimuhat0,dimuhat0),zeros(dimuhat0,dimuhat0),G];
gamma_mean_trace_local=zeros(dimuhat,N);
covph=zeros(3*dimuhat,3*dimuhat,N/gap);
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
 uhatmesh=[invG1;invG2;invG3]*uhatnew;%%%%% transform to phy domain
 if max(abs(imag( uhatmesh(:))))>10^(-6)
     disp('error of transform to phy domain')
 else
     uhatmesh=real(uhatmesh);
 end
   [avexens,xensnew]=compute_ensave_bd(xensnew,bdryh,MM);
  [aveyens,yensnew]=compute_ensave_bd(yensnew,bdryh,MM);
 temp_obs=[xobs(:,i);yobs(:,i)]; 
idbdobsx=find(abs(avexens-temp_obs(1:L))>2*pi-2*bdryh);
idbdobsy=find(abs(aveyens-temp_obs(1+L:end))>2*pi-2*bdryh);
temp_obs(idbdobsx)=temp_obs(idbdobsx)-2*pi;
temp_obs(idbdobsy+L)=temp_obs(idbdobsy+L)-2*pi;

   u_prior=[xensnew;yensnew;uhatmesh];
    u_mean_prior=sum(u_prior,2)/MM;
    U=u_prior - u_mean_prior * ones(1, MM);
    U=sqrt(1+r)*U;
    V=U(1:2*L,:);%%%% observation residual 
aveuhat=sum(uhatmesh,2)/MM;  

      Ux=U(1:L,:);Uy=U(1+L:2*L,:);
    Uhat=U(1+2*L:end,:);
  %% loop for drift
driftid=get_local_driftID(gamma,avexens,aveyens,temp_obs(1:L),temp_obs(1+L:end));%%%% local id for all drift
localaveVx=cell(L,1);localaveVy=cell(L,1);%% \barx^f_l
localVx=cell(L,1);localVy=cell(L,1);%%%% X_l^f
xobslocal=cell(L,1);yobslocal=cell(L,1);%%% y^o_l
for idrift=1:L
    localindx=driftid{idrift};
localaveVx{idrift}=avexens(localindx);localVx{idrift}=V(localindx,:);
localaveVy{idrift}=aveyens(localindx);localVy{idrift}=V(localindx+L,:);
xobslocal{idrift}=temp_obs(localindx);yobslocal{idrift}=temp_obs(localindx+L);
end

ensxup=zeros(L,MM);ensyup=zeros(L,MM);
avexup=zeros(L,1);aveyup=zeros(L,1);
for idrift=1:L
% for idrift=1:L 
    if  size(driftid{idrift},1)>0
    idrift1=idrift+L;
    Ylf=[localVx{idrift};localVy{idrift}];
    aveylf=[localaveVx{idrift};localaveVy{idrift}];
    yol=[xobslocal{idrift};yobslocal{idrift}];
    Q=(MM-1)*speye(MM)/(1+r)+1/(sigma_obs.^2)*(Ylf'*Ylf);
%     norm(Q)
%     Pla=pinv(Q);
        [V1, D] = eig(Q/2+Q'/2); 
%     [V1, D] = eig(Q); 
    Pla=V1*diag(1./diag(D))*V1';%%%a=Pla-pinv(Q);max(abs(a(:)))      
%     Wlatemp=sqrt(MM-1)*V1*diag(sqrt(1./diag(D)));
    Wlatemp=sqrt(MM-1)*V1*diag(sqrt(1./diag(D)))*V1';%%%%T  
    avewla=1/(sigma_obs.^2)*Pla*Ylf'*(yol-aveylf);
    Wla=Wlatemp+avewla*ones(1,MM);
    ensxup(idrift,:)=avexens(idrift)*ones(1,MM)+Ux(idrift,:)*Wla;
    ensyup(idrift,:)=aveyens(idrift)*ones(1,MM)+Uy(idrift,:)*Wla;
    avexup(idrift)=avexens(idrift)+Ux(idrift,:)*avewla;
    aveyup(idrift)=aveyens(idrift)+Uy(idrift,:)*avewla;
  else
       ensxup(idrift,:)=xensnew(idrift,:); ensyup(idrift,:)=yensnew(idrift,:);
       avexup(idrift)=avexens(idrift);aveyup(idrift)=aveyens(idrift);
    end
end
xens=ensxup;
yens=ensyup;%%%% update for next predict
%% loop for uhat

idmeshdrift=get_local_vdID((2*K_max+1),gamma,temp_obs(1:L),temp_obs(1+L:end),2*pi,0);%%%% local id for all meshgrid
localaveVx=cell(dimuhat,1);localaveVy=cell(dimuhat,1);%% \bary^f_l
localVx=cell(dimuhat,1);localVy=cell(dimuhat,1);%%%% Y_l^f
xobslocal=cell(dimuhat,1);yobslocal=cell(dimuhat,1);%%% y^o_l
for imesh=1:dimuhat0
    localindx=idmeshdrift{imesh};
     allnobs(imesh,i/gap)=length(localindx);
localaveVx{imesh}=avexens(localindx);localVx{imesh}=V(localindx,:);
localaveVy{imesh}=aveyens(localindx);localVy{imesh}=V(localindx+L,:);
xobslocal{imesh}=temp_obs(localindx);yobslocal{imesh}=temp_obs(localindx+L);
end
ensuhatup=zeros(dimuhat,MM);
aveuhatup=zeros(dimuhat,1);
% pause()

for imesh=1:dimuhat0
% parfor imesh=1:dimuhat0
    if  size(idmeshdrift{imesh},2)>0
      Ylf=[localVx{imesh};localVy{imesh}];
    aveylf=[localaveVx{imesh};localaveVy{imesh}];
    yol=[xobslocal{imesh};yobslocal{imesh}];
    Q=(MM-1)*speye(MM)/(1+r)+1/(sigma_obs.^2)*(Ylf'*Ylf);
%     Pla=pinv(Q);
        [V1, D] = eig(Q/2+Q'/2); 
%     [V1, D] = eig(Q); 
    Pla=V1*diag(1./diag(D))*V1';%%%a=Pla-pinv(Q);max(abs(a(:)))  
% %     Wlatemp=sqrt(MM-1)*V1*diag(sqrt(1./diag(D)));
    Wlatemp=sqrt(MM-1)*V1*diag(sqrt(1./diag(D)))*V1';
    avewla=1/(sigma_obs.^2)*Pla*Ylf'*(yol-aveylf);
% % % s=V*diag(sqrt(1./diag(D)));a=s*s'-Pla;max(abs(a(:)))
    Wla=Wlatemp+avewla*ones(1,MM);  
ensuhatup([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0],:)=aveuhat([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0])*ones(1,MM)+...
        Uhat([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0],:)*Wla;
aveuhatup([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0])=aveuhat([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0])+...
    Uhat([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0],:)*avewla;

    else
       ensuhatup([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0],:)=uhatmesh([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0],:);
       aveuhatup([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0])=aveuhat([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0]);
    end
end %% end of analysi 

  Gu=G*ensuhatup(1:dimuhat0,:);%%%% temp transform to F domain
  Gv=G*ensuhatup(1+dimuhat0:2*dimuhat0,:);
% % %    Gu=G*uhatmesh(1:dimuhat0,:);%%%% temp transform to F domain
% % %   Gv=G*uhatmesh(1+dimuhat0:2*dimuhat0,:);
 
  u_mean_post=sum(ensuhatup,2)/MM;
    Upost=ensuhatup - u_mean_post * ones(1, MM);
covph(:,:,i/gap)=Upost*Upost';

% uhat=zeros(dimuhat,MM);
uhat=uhatnew;
for id=2:dimuhat0
%     if length(idmeshdrift{id})==0%%%% of no observation then no update
%     disp('no obs')
%     uhat(id,:)=uhatnew(id,:);
%     else
    if rk(1,id)==0
      uhat(id,:)=Gv(id,:)/rk(2,id);  
    else
     uhat(id,:)=Gu(id,:)/rk(1,id);   
    end
%     end
end

    uhat(1,:)=0;
%     pause()
 gamma_mean_trace_local(:,i)=mean(uhat,2);  
%     blowuptest(real(gamma_mean0),10^6);
if max(real(uhat(:)))>10^5
    disp('blow up,return')
    i
    return
end
else%%% no update
     uhat=uhatnew;
     xens= xensnew;
     yens= yensnew;
    xens=mod(xens,2*pi);
    yens=mod(yens,2*pi);
end
end
fprintf('gap is %d\n',gap);
fprintf('dt gap is %2.3f\n',gap*dt);
fprintf('mean observation is %2.2f\n',mean(allnobs(:)));
% % figure();plot(dt:dt:T,real(u_hat(ind(1),:)),'r','linewidth',1.5);hold on;
% % plot(dt:dt:T,real(gamma_mean_trace_local(ind(1),:)),'b','linewidth',1.5);
% % setgca(16);title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
% % xlabel('t','fontsize',16)
rmscc(real(gamma_mean_trace_local(ind(1),gap:gap:end)),real(u_hat(ind(1),gap:gap:end)),1);

