%%%%% localized conditional gaussian
close all;
updatecov=1;
gamma=28;fprintf('gamma is %2.2f\n',gamma);
disp('localization GB without r');
F=fft(eye(2*K_max+1));
fft2mat = kron(F,F);invfft2mat=fft2mat\eye((2*K_max+1)^2);
fixedcov=diag((sigmak2exact.^2)./(dkexact+sqrt(dkexact.^2+L/sig_ex/sig_ex.*sigmak2exact.^2)));
fixedcov(1,1)=1;
gamma_cov0=fixedcov;
% gamma_cov0 = eye(Dim_Yr)*0.01; 
gamma_mean0=zeros(Dim_Yr,1);
gamma_mean_trace_local=zeros(Dim_Y,N);%%% store all
gamma_mean_trace_local(redind,1) = gamma_mean0;
disp('localized data assimilation (filter)......')
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
    
%   idmeshdrift=get_local_vdID((2*K_max+1),gamma,xexact(:,i-1),yexact(:,i-1),L1,0);  
  idmeshdrift=get_local_vdID_L((2*K_max+1),6,xexact(:,i),yexact(:,i),L1);

  gamma_mean=zeros(length(idmeshdrift),1);
  covph=zeros(length(idmeshdrift),1);
    for j=1:length(idmeshdrift)
%         idx=1:2*L;
idx=[idmeshdrift{j},idmeshdrift{j}+L];
idx=randperm(L, 1);idx=[idx,idx+L];
% idx=[j,j+L];
            Lo=invfft2mat(j,:);
    gamma_mean(j) = Lo*(gamma_mean0 + (a0 + a1 * gamma_mean0) * dt) + ...
        (Lo*gamma_cov0 * A1(idx,:)') * (invBoBdiag(idx).*  (diff_xy(idx) - A0(idx)*dt-A1(idx,:) * gamma_mean0 * dt));
        if imag(abs(gamma_mean(j)))>10^(-6)
        disp('error in localization (complex)')
        else
        gamma_mean(j)=real( gamma_mean(j));
        end
covph(j) = Lo*(gamma_cov0 + (a1 * gamma_cov0 + ...
gamma_cov0 * a1' + b1b1t)*dt)*Lo' - ((Lo*gamma_cov0 * A1(idx,:)') * invBoB(idx,idx) * (Lo*gamma_cov0 * A1(idx,:)')') * dt;     
    end
    gamma_mean0=fft2(reshape(gamma_mean,2*K_max+1,2*K_max+1));gamma_mean0=gamma_mean0(:);
% %     gamma_meant=fft2mat*gamma_mean;

 %%%%% physical  domain global update
% % covph1 = invfft2mat*(gamma_cov0 + (a1 * gamma_cov0 + ...
% % gamma_cov0 * a1' + b1b1t)*dt)*invfft2mat' - ...
% % ((invfft2mat*gamma_cov0 * A1') * invBoB * (invfft2mat*gamma_cov0 * A1')') * dt;     

% gamma_cov0 = (gamma_cov0 + (a1 * gamma_cov0 + ...
% gamma_cov0 * a1' + b1b1t)*dt) - ((gamma_cov0 * A1') * invBoB * ...
% (gamma_cov0 * A1')') * dt;     %%%%% fourier domain update
% 

gamma_mean_trace_local(redind,i) = gamma_mean0;
if updatecov
gamma_cov0=fft2mat*diag(covph)*fft2mat';%gamma_cov0=real(gamma_cov0);%%%physical  domain local update
gamma_cov0(1,1)=1;
end
% gamma_cov0=fft2mat*(covph1)*fft2mat';gamma_cov0=real(gamma_cov0);%%% physical  domain global update

end
 [rmstimeloc,rmsfilloc,ccfilloc ]=compute_velerrorGBsw(u_hat(:,gap:gap:end),gamma_mean_trace_local(:,gap:gap:end),rk,kk,L1,Dim_Grid,3);
 errfillocal=[mean(rmstimeloc(3:end));mean(rmsfilloc);mean(ccfilloc)];

% figure();plot(dt:dt:T,real(u_hat(ind,:)),'r','linewidth',1.5);hold on;
% plot(dt:dt:T,real(gamma_mean_trace_local(ind,:)),'b','linewidth',1.5);
% setgca(16);title(['Filter, mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)
% xlabel('t','fontsize',16)
rmscc(real(gamma_mean_trace(ind(1),:)),real(u_hat(ind(1),:)),1);
disp('......')
rmscc(real(gamma_mean_trace_local(ind(1),:)),real(u_hat(ind(1),:)),1);

% rmscc(real(gamma_mean_trace(ind(1),:)),real(gamma_mean_trace_local(ind(1),:)),1);

return
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