%%%%% localized conditional gaussian
%%%%% use etkf to get the covariance matrix
close all;
updatecov=0;dimuhat0=dimuhat;localupR=1;MM=200;
fprintf('ensemble size is %d\n',MM);r=0;sigma_obs=0.1;
gamma=28;fprintf('gamma is %2.2f\n',gamma);
% F=fft(eye(2*K_max+1));
% fft2mat = kron(F,F);invfft2mat=fft2mat\eye((2*K_max+1)^2);
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
fixedcov=diag((sigmak2exact.^2)./(dkexact+sqrt(dkexact.^2+L/sig_ex/sig_ex.*sigmak2exact.^2)));
fixedcov(1,1)=1;
gamma_cov0=fixedcov;
% gamma_cov0 = eye(Dim_Yr)*0.01; 
gamma_mean0=zeros(Dim_Yr,1);
gamma_mean_trace_local=zeros(Dim_Y,N);%%% store all
gamma_mean_trace_local(redind,1) = gamma_mean0;
allnobs=zeros(dimuhat0,N);
Ro=(sigma_obs).^2*speye(2*L);
H=speye(2*L,2*L+dimuhat*3);nfloe=L;
pricov=zeros(dimuhat*3,dimuhat*3,N);
postcov=zeros(dimuhat*3,dimuhat*3,N);

disp('localized data assimilation (filter)......')
for i = 3:N
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
  u=zeros(length(idmeshdrift),1);
  v=zeros(length(idmeshdrift),1);
  covph=zeros(length(idmeshdrift)*3,1);
    for j=1:length(idmeshdrift)
        allnobs(j,i)=length(idmeshdrift{j});
%         idx=1:2*L;
% if length(idmeshdrift{j})==0
% %     disp('no obs')
% else
%     idmeshdrift{j}=1:L;
% end
idx=[idmeshdrift{j},idmeshdrift{j}+L];
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
   covph(j+2*dimuhat0) = Leta*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Leta' - ((Leta*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Leta*gamma_cov0 * A1(idx,:)')') * dt;     

        if imag(abs(v(j)))>10^(-6) || imag(abs(u(j)))>10^(-6)
        disp('error in localization (complex)')
        else
        v(j)=real( v(j));u(j)=real( u(j));
        end
    end
    %% transform to F domain
        %%% transform to Fourier domain
Gu=G*u;
Gv=G*v;
etkf_analysis_gb
for id=2:dimuhat0
    if length(idmeshdrift{id})==0%%%% of no observation then no update
% %     disp('no obs')
    else
    ka=kk0(1,id);kb=kk0(2,id);
    if rk(1,id)==0
      gamma_mean0(id)=Gv(id)/rk(2,id);  
    else
       gamma_mean0(id)=Gu(id)/rk(1,id);   
    end
    end
end
gamma_mean0(1)=0;
 %%%%% physical  domain global update
% % covph1 = invfft2mat*(gamma_cov0 + (a1 * gamma_cov0 + ...
% % gamma_cov0 * a1' + b1b1t)*dt)*invfft2mat' - ...
% % ((invfft2mat*gamma_cov0 * A1') * invBoB * (invfft2mat*gamma_cov0 * A1')') * dt;     
% 
% 
if max(real(gamma_mean0))>10^5
    disp('blow up,return')
    i
    return
end
gamma_mean_trace_local(redind,i) = gamma_mean0;
if localupR
gamma_cov0=Gm*diag(covph)*Gm';upcov=1;
% gamma_cov0=Gm*covph1*Gm';upcov=1;
else
% gamma_cov0 = (gamma_cov0 + (a1 * gamma_cov0 + ...
% gamma_cov0 * a1' + b1b1t)*dt) - ((gamma_cov0 * A1') * invBoB * ...
% (gamma_cov0 * A1')') * dt;
    upcov=2;     %%%%% fourier global domain update
end
% gamma_cov0=real(gamma_cov0);%%%physical  domain local update
% end
% % gamma_cov0=fft2mat*(covph1)*fft2mat';gamma_cov0=real(gamma_cov0);%%% physical  domain global update

end
if upcov==1
    disp('update cov locally in physical domain')
else
    disp('update cov global in Fourier domain')
end
fprintf('gamma is %2.2f\n',gamma);
fprintf('mean observation is %2.2f\n',mean(s(:)));
fprintf('gamma is %2.2f\n',gamma);

figure();plot(dt:dt:T,real(u_hat(ind,:)),'r','linewidth',1.5);hold on;
plot(dt:dt:T,real(gamma_mean_trace_local(ind,:)),'b','linewidth',1.5);
setgca(16);title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
xlabel('t','fontsize',16)
rmscc(real(gamma_mean_trace(ind(1),3:end)),real(u_hat(ind(1),3:end)),1);
disp('......')
rmscc(real(gamma_mean_trace_local(ind(1),3:end)),real(u_hat(ind(1),3:end)),1);

% rmscc(real(gamma_mean_trace_local(ind(1),:)),real(gamma_mean_trace(ind(1),:)),1);

return
gap=200;ix=2;iy=2;
load dataex;load data25;load data30;load datainf;load data21
figure();plot(dt:dt*gap:T,real(dataex(gap:gap:end)),'b','linewidth',1.5);hold on;
plot(dt:dt*gap:T,real(datainf(gap:gap:end)),'r','linewidth',1.5);hold on
plot(dt:dt*gap:T,real(data30(gap:gap:end)),'k','linewidth',1.5);hold on
plot(dt:dt*gap:T,real(data25(gap:gap:end)),'g','linewidth',1.5);
plot(dt:dt*gap:T,real(data21(gap:gap:end)),'c','linewidth',1.5);
setgca(16);
lgnd=legend('Exact','No localization','r=3','r=2.5','r=2.1');
set(lgnd,'FontSize',15);

title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',20)
xlabel('t','fontsize',20)




gap=200;ix=2;iy=2;
load dataex;load data25;load data30;load datainf;load data35
figure();plot(dt:dt*gap:T,real(dataex(gap:gap:end)),'b','linewidth',1.5);hold on;
plot(dt:dt*gap:T,real(datainf(gap:gap:end)),'r','linewidth',1.5);hold on
plot(dt:dt*gap:T,real(data35(gap:gap:end)),'k','linewidth',1.5);hold on
plot(dt:dt*gap:T,real(data30(gap:gap:end)),'g','linewidth',1.5);
plot(dt:dt*gap:T,real(data25(gap:gap:end)),'c','linewidth',1.5);
setgca(16);
lgnd=legend('Exact','No localization','r=3.5','r=3.0','r=2.5');
set(lgnd,'FontSize',15);

title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',20)
xlabel('t','fontsize',20)


%%%%
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