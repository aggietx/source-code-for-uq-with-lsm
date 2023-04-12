%%%%% localized conditional gaussian
%%%% shallow water%% no update for no obs
gamma=30;localupR=1;fprintf('gamma is %2.2f\n',gamma);
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
% Gm=repmat(G,3,3);
sigmak2exact=[sigmak2exactb;sigmak2exactg;sigmak2exactg];
fixedcov=diag((sigmak2exact.^2)./(dkexact+sqrt(dkexact.^2+L/sig_ex/sig_ex.*sigmak2exact.^2)));
fixedcov(1,1)=1;fixedcov(1+dimuhat0,1+dimuhat0)=1;fixedcov(1+dimuhat0*2,1+dimuhat0*2)=1;
gamma_cov0=fixedcov;
gamma_mean0=zeros(Dim_Yr,1);
gamma_mean_trace_local=zeros(Dim_Y,N);%%% store all
gamma_mean_trace_local(redind,1) = gamma_mean0;

allnobs=zeros(dimuhat0,N);
gap=100;
Dt=gap*dt;Np=T/Dt;
allcc=zeros(Np,1);allrms=zeros(Np,1);

disp('localized data assimilation (filter)......')
rng(10)
for i = 2:N
    if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
    end
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
    
    A1=[G1;G2];%A1=A1(:,redind);
    
  idmeshdrift=get_local_vdID((2*K_max+1),gamma,xexact(:,i),yexact(:,i),L1,0); 
%   idmeshdrift=get_local_vdID_L((2*K_max+1),6,xexact(:,i),yexact(:,i),L1);

    u=zeros(length(idmeshdrift),1);
  v=zeros(length(idmeshdrift),1);
  eta=zeros(length(idmeshdrift),1);
covph=zeros(length(idmeshdrift)*3,1);  
    for j=1:length(idmeshdrift)
%         allnobs(j,i)=length(idmeshdrift{j});
%         idx=1:2*L;
% idx=[idmeshdrift{j},idmeshdrift{j}+L];
idx=randperm(L, 12);idx=[idx,idx+L];
% idx=1:2*L;
% if isempty(idmeshdrift{j})
% %     disp('no obs')
% else
%     idmeshdrift{j}=1:L;
% end
 allnobs(j,i)=length(idx)/2;
            Lu=invG1(j,:);
   u(j) = Lu*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lu*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
   covph(j) = Lu*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lu' - ((Lu*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lu*gamma_cov0 * A1(idx,:)')') * dt;     
            
    Lv=invG2(j,:);
    v(j) = Lv*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lv*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
   covph(j+dimuhat0) = Lv*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lv' - ((Lv*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lv*gamma_cov0 * A1(idx,:)')') * dt;     
    
               Leta=invG3(j,:);
   eta(j) = Leta*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Leta*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
   covph(j+2*dimuhat0) = Leta*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Leta' - ((Leta*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Leta*gamma_cov0 * A1(idx,:)')') * dt;     

    if imag(abs(u(j)))>10^(-6) || imag(abs(v(j)))>10^(-6) || imag(abs(eta(j)))>10^(-6)
        disp('error in localization (complex)')
        else
       u(j)=real( u(j));eta(j)=real( eta(j));v(j)=real( v(j));
    end
% %     gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + ...
% % gamma_cov0 * a1' + b1b1t - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
    end
  ub=G*u;
  vb=G*v;
  etab=G*eta;
  
%   gamma_mean0=zeros(dimuhat,1);
    for ii=2:length(kk0)
        if length(idmeshdrift{ii})==0%%%% of no observation then no update
% %     disp('no obs')            
        else
       ind=[ii,ii+dimuhat0,ii+2*dimuhat0];
       ux=ub(ii);vx=vb(ii);etax=etab(ii);
       rmatrix=[rk(:,ind);rketa(:,ind)];
       gamma_mean0(ind)=rmatrix\[ux;vx;etax];
%         end
        end
    end  
    gamma_mean0([1,1+dimuhat0,1+2*dimuhat0])=0;
%     blowuptest(real(gamma_mean0),10^6);
if max(real(gamma_mean0))>10^5
    disp('blow up,return')
    i
    return
end

gamma_mean_trace_local(redind,i) = gamma_mean0;

if localupR
% gamma_cov0=Gm*diag(covph)*Gm';upcov=1;
gamma_cov0temp=Gm*diag(covph)*Gm';upcov=1;
gamma_cov0=zeros(size(gamma_cov0temp));
    for ii=2:length(kk0)
%         if length(idmeshdrift{ii})==0%%%% of no observation then no update
% %     disp('no obs')            
%         else
       ind=[ii,ii+dimuhat0,ii+2*dimuhat0];
       rmatrix=[rk(:,ind);rketa(:,ind)];
      gamma_cov0(ind,ind)=inv(rmatrix)*gamma_cov0temp(ind,ind)*inv(rmatrix');  
%     gamma_cov0(ind,ind)=rmatrix*gamma_cov0temp(ind,ind)*rmatrix';            
% gamma_cov0(ind,ind)=gamma_cov0temp(ind,ind); 
%         end
    end
    gamma_cov0([1,1+dimuhat0,1+2*dimuhat0],[1,1+dimuhat0,1+2*dimuhat0])=eye(3);
else
   gamma_cov0 = (gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt) - ((gamma_cov0 * A1') * invBoB * ...
(gamma_cov0 * A1')') * dt;  upcov=2;   %%%%% fourier global domain update
 
end
gamma_cov0=(gamma_cov0+gamma_cov0.')/2;
% % gamma_cov0=real(gamma_cov0);%%%physical  domain local update
% % end
if mod(i,gap)==0
    [ux,uy]=computeGBvel(rk,kk,L1,Dim_Grid,u_hat(:,i));
    [uxest,uyest]=computeGBvel(rk,kk,L1,Dim_Grid,gamma_mean0);
[rmsvel,ccvel]=rmscc(sqrt(uxest.^2+uyest.^2),sqrt(ux.^2+uy.^2),0);
allcc(i/gap)=ccvel;
allrms(i/gap)=rmsvel;
end
end
if upcov==1
    disp('update cov locally in physical domain')
else
    disp('update cov global in Fourier domain')
end
fprintf('lcg, mean rms, cc are  %2.2f  %2.2f\n',mean(allrms),mean(allcc));
% save allcc allcc;
% save allrms allrms
s=allnobs(:,4:end);
fprintf('mean observation is %2.2f\n',mean(s(:)));

fprintf('gamma is %2.2f\n',gamma);
% figure();plot(dt:dt:T,real(u_hat(ind(1),:)),'r','linewidth',1.5);hold on;
% plot(dt:dt:T,real(gamma_mean_trace_local(ind(1),:)),'b','linewidth',1.5);
% setgca(16);title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
% xlabel('t','fontsize',16)
rmscc(real(gamma_mean_trace(ind(1),:)),real(u_hat(ind(1),:)),1);
rmscc(real(gamma_mean_trace(ind(2),:)),real(u_hat(ind(2),:)),1);

disp('....')
rmscc(real(gamma_mean_trace_local(ind(1),:)),real(u_hat(ind(1),:)),1);
rmscc(real(gamma_mean_trace_local(ind(2),:)),real(u_hat(ind(2),:)),1);


return
gap=50;ix=2;iy=2;T=20;
load dataex;load data25;load data35;load datainf;load data350
figure();subplot(2,1,1)
plot(dt:dt*gap:T,real(dataex(1,gap:gap:end)),'b','linewidth',1.5);hold on;
plot(dt:dt*gap:T,real(datainf(1,gap:gap:end)),'r','linewidth',1.5);hold on
plot(dt:dt*gap:T,real(data25(1,gap:gap:end)),'k','linewidth',1.5);hold on
plot(dt:dt*gap:T,real(data35(1,gap:gap:end)),'g','linewidth',1.5);
plot(dt:dt*gap:T,real(data350(1,gap:gap:end)),'c','linewidth',1.5);
setgca(16);
% lgnd=legend('Exact','No localization','r=3','r=2.5','r=2.1');
% set(lgnd,'FontSize',15);
title(['Filter, GB mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',20)
% xlabel('t','fontsize',20)
subplot(2,1,2)
plot(dt:dt*gap:T,real(dataex(2,gap:gap:end)),'b','linewidth',1.5);hold on;
plot(dt:dt*gap:T,real(datainf(2,gap:gap:end)),'r','linewidth',1.5);hold on
plot(dt:dt*gap:T,real(data25(2,gap:gap:end)),'k','linewidth',1.5);hold on
plot(dt:dt*gap:T,real(data35(2,gap:gap:end)),'g','linewidth',1.5);
plot(dt:dt*gap:T,real(data350(2,gap:gap:end)),'c','linewidth',1.5);
setgca(16);
lgnd=legend('Exact','No localization (no transform) ','r=2.5','r=3.5','r=35 (transform)');
set(lgnd,'FontSize',15);
title(['Filter, gr mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',20)
xlabel('t','fontsize',20)
return

dimt=20;
F=fft(eye(2*dimt+1));
fft2mat = kron(F,F);invfft2mat=fft2mat\eye((2*dimt+1)^2);

a=randn((2*dimt+1)^2,1);

tic;aa=reshape(fft2mat*a,2*dimt+1,2*dimt+1);toc
tic;bb=fft2(reshape(a,2*dimt+1,2*dimt+1));toc
 norm(aa-bb)
 
 tic;dd=reshape(invfft2mat*aa(:),2*dimt+1,2*dimt+1);toc
 tic;cc=ifft2(reshape(bb,2*dimt+1,2*dimt+1));toc
  norm(dd-cc)