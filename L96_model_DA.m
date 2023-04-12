clear;close all;
I=40;J=5;Linewid=2;fonsize=18;
T = 150;
dt = 2e-2;
N = round(T/dt);
xlim=50:dt:150;Tacf=15;
% range=floor(50/dt):floor(150/dt);
range=1:N;
u_truth=zeros(I,N);
v_truth=zeros(I*J,N);%%% (i-1)*j+1:i*j

d_i_bar=ones(1,I)+0.7*cos(2*pi*(1:I)/J);
gamma_i=0.1*ones(1,I)+0.25*cos(2*pi*(1:I)/J);
gamma_0=0.1*ones(1,J)+0.25*cos(2*pi*(1:J)/J);
% gamma_j=repmat(gamma_0,I,1);gamma_j=gamma_j(:);
gamma_j=repmat(gamma_i',J,1);
sigma_u=1;
sigma_ij=zeros(I*J,1);sigma_ij((1-1)*I+1:1*I)=0.5;sigma_ij((2-1)*I+1:2*I)=0.2;
sigma_ij((3-1)*I+1:3*I)=0.1;sigma_ij((4-1)*I+1:4*I)=0.1;sigma_ij((5-1)*I+1:5*I)=0.1;

d_vij=zeros(I*J,1);d_vij((1-1)*I+1:1*I)=0.2;d_vij((2-1)*I+1:2*I)=0.5;
d_vij((3-1)*I+1:3*I)=1;d_vij((4-1)*I+1:4*I)=2;d_vij((5-1)*I+1:5*I)=5;
% d_vij=5*ones(I*J,1);
F=8;

xlimset=[50,150];
for i = 2:N
   v_truthtemp=v_truth(:,i-1);v_truthtemp=reshape(v_truthtemp,I,J);
   u_truth(1,i)= u_truth(1,i-1)+dt*(u_truth(I,i-1).*(u_truth(2,i-1)-u_truth(I-1,i-1))-d_i_bar(1)'.*u_truth(1,i-1)+F+...
       u_truth(1,i-1).*gamma_i(1)'.*sum(v_truthtemp(1,:))) + sigma_u * sqrt(dt) * randn;  
   u_truth(2,i)= u_truth(2,i-1)+dt*(u_truth(1,i-1).*(u_truth(3,i-1)-u_truth(I,i-1))-d_i_bar(2)'.*u_truth(2,i-1)+F+...
       u_truth(2,i-1).*gamma_i(2)'.*sum(v_truthtemp(2,:)))+ sigma_u * sqrt(dt) * randn;   
   u_truth(3:I-1,i)=u_truth(3:I-1,i-1)+dt*(u_truth(2:I-2,i-1).*(u_truth(4:I,i-1)-u_truth(1:I-3,i-1))-d_i_bar(3:I-1)'.*u_truth(3:I-1,i-1)+F+...
       u_truth(3:I-1,i-1).*gamma_i(3:I-1)'.*sum(v_truthtemp(3:I-1,:),2) )+sigma_u * sqrt(dt) * randn(I-3,1);  
   u_truth(I,i)= u_truth(I,i-1)+dt*(u_truth(I-1,i-1).*(u_truth(1,i-1)-u_truth(I-2,i-1))-d_i_bar(I)'.*u_truth(I,i-1)+F+...
       u_truth(I,i-1).*gamma_i(I)'.*sum(v_truthtemp(I,:)))+ sigma_u * sqrt(dt) * randn;  
   
   v_truth(:,i)=v_truth(:,i-1)+dt*(-d_vij.*v_truth(:,i-1) -gamma_j.*repmat(u_truth(:,i-1).^2,J,1))+ sqrt(dt) *sigma_ij.* randn(I*J,1);
end

a=u_truth(:,1:floor(30/dt))';%a=a(end:-1:1,:);
% figure();imagesc(a);colorbar
% figure();plot(xlim,v_truth(21,range));%xlim([50,150])
% return
disp('data assi filter..');
gamma_mean0 = zeros(I*J,1);
gamma_cov0 = zeros(I*J,I*J);
gamma_mean_trace=zeros(I*J,N);
gamma_cov_trace=zeros(I^2*J^2,N);
for i = 2:N
    u0 = u_truth(:,i-1);
    u = u_truth(:,i);
    a0=-gamma_j.*repmat(u0.^2,J,1);
    a1=sparse(1:I*J,1:I*J,-d_vij);
    b1=sparse(1:I*J,1:I*J,sigma_ij);
    invBoB=sparse(1:I,1:I,1./(sigma_u^2)*ones(I,1));
    A0=[u0(I).*(u0(2)-u0(I-1))-d_i_bar(1)'.*u0(1)+F;u0(1).*(u0(3)-u0(I))-d_i_bar(2)'.*u0(2)+F;...
        u0(2:I-2).*(u0(4:I)-u0(1:I-3))-d_i_bar(3:I-1)'.*u0(3:I-1)+F;u0(I-1).*(u0(1)-u0(I-2))-d_i_bar(I)'.*u0(I)+F];
    A1=zeros(I,I*J);
    for ii=1:I
        A1(ii,ii:I:end)=u0(ii).*gamma_i(ii)*ones(J,1);
    end
    gamma_mean = gamma_mean0 + (a0 + a1* gamma_mean0) * dt + ...
    (gamma_cov0 * A1') * invBoB * (u-u0 - A0*dt-A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' ...
        - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt; 
    
    gamma_mean_trace(:,i) = gamma_mean;
    gamma_cov_trace(:,i) = reshape(gamma_cov,I^2*J^2,[]);

    gamma_mean0 = gamma_mean;
    gamma_cov0 = gamma_cov;
    
end
% ind=21;
% plot2line(v_truth(ind,range),gamma_mean_trace(ind,range))
% ind=21+40;
% plot2line(v_truth(ind,range),gamma_mean_trace(ind,range))
% ind=21+80;
% plot2line(v_truth(ind,range),gamma_mean_trace(ind,range))
disp('data assi smoother..');
%% backward
gamma_mean0_smoother = gamma_mean0;
gamma_cov0_smoother = gamma_cov0 ;
gamma_mean_trace_smoother=zeros(I*J,N);
gamma_mean_trace_smoother(:,end)=gamma_mean0;

for i = N:-1:2
    u0 = u_truth(:,i);
 R=gamma_cov_trace(:,i);R=reshape(R,I*J,I*J);invR=R\speye(I*J);
    a0=-gamma_j.*repmat(u0.^2,J,1);
    a1=sparse(1:I*J,1:I*J,-d_vij);
    b1=sparse(1:I*J,1:I*J,sigma_ij);

        gamma_mean_smoother=gamma_mean0_smoother+dt*(-a0-a1*gamma_mean0_smoother+...
            (b1 * b1')*invR*(gamma_mean_trace(:,i)-gamma_mean0_smoother));
    gamma_cov_smoother=gamma_cov0_smoother-dt*( (a1+(b1*b1')*invR)*gamma_cov0_smoother+...
        gamma_cov0_smoother*(a1'+b1*b1'*invR)-b1*b1' );
    
        gamma_mean_trace_smoother(:,i-1) = gamma_mean_smoother;

    gamma_mean0_smoother = gamma_mean_smoother;
    gamma_cov0_smoother = gamma_cov_smoother;

    
end
% ind=21;
% plot2line(v_truth(ind,range),gamma_mean_trace_smoother(ind,range))
% ind=21+40;
% plot2line(v_truth(ind,range),gamma_mean_trace_smoother(ind,range))
% ind=21+80;
% plot2line(v_truth(ind,range),gamma_mean_trace_smoother(ind,range))

%% backward sample
gamma_mean0_smoother_sample = gamma_mean0;
gamma_mean_trace_smoother_sample=zeros(I*J,N);
gamma_mean_trace_smoother_sample(:,end)=gamma_mean0;

for i = N:-1:2
    u0 = u_truth(:,i);
 R=gamma_cov_trace(:,i);R=reshape(R,I*J,I*J);invR=R\speye(I*J);
    a0=-gamma_j.*repmat(u0.^2,J,1);
    a1=sparse(1:I*J,1:I*J,-d_vij);
    b1=sparse(1:I*J,1:I*J,sigma_ij);

        gamma_mean_smoother_sample=gamma_mean0_smoother_sample+dt*(-a0-a1*gamma_mean0_smoother_sample+...
            (b1 * b1')*invR*(gamma_mean_trace(:,i)-gamma_mean0_smoother_sample))+sqrt(dt)*b1*randn(I*J,1);
        
        gamma_mean_trace_smoother_sample(:,i-1) = gamma_mean_smoother_sample;

    gamma_mean0_smoother_sample = gamma_mean_smoother_sample;

    
end
return
ind=21;
% ind=21+40;
% ind=21+80;
% plot2line(v_truth(ind,range),gamma_mean_trace_smoother_sample(ind,range))
% range=1:N;
figure();[~,r1]=ksdensity(u_truth(ind,:));plot(r1,ksdensity(u_truth(ind,range)),'r','Linewidth',2);
plotpdf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),gamma_mean_trace_smoother_sample(ind,:),range)
plotacf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),...
    gamma_mean_trace_smoother_sample(ind,:),1:N,Tacf,dt)
% rnorm(gamma_mean_trace(ind,ranget),v_truth(ind,ranget));
% rnorm(gamma_mean_trace_smoother(ind,ranget),v_truth(ind,ranget));
% rnorm(gamma_mean_trace_smoother_sample(ind,ranget),v_truth(ind,ranget));
% return
ind=21+40;
plotpdf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),gamma_mean_trace_smoother_sample(ind,:),range)
plotacf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),...
    gamma_mean_trace_smoother_sample(ind,:),1:N,Tacf,dt)

ind=21+80;
plotpdf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),gamma_mean_trace_smoother_sample(ind,:),range)
plotacf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),...
    gamma_mean_trace_smoother_sample(ind,:),1:N,Tacf,dt)

ind=21+120;
plotpdf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),gamma_mean_trace_smoother_sample(ind,:),range)
plotacf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),...
    gamma_mean_trace_smoother_sample(ind,:),1:N,Tacf,dt)

ind=21+160;
plotpdf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),gamma_mean_trace_smoother_sample(ind,:),range)
plotacf(v_truth(ind,:),gamma_mean_trace(ind,:),gamma_mean_trace_smoother(ind,:),...
    gamma_mean_trace_smoother_sample(ind,:),1:N,Tacf,dt)

return
e=zeros(I*J,1);
for ind=1:I*J
    close all
    ind
  plot2line(v_truth(ind,:),gamma_mean_trace(ind,:))
  e(ind)=norm((v_truth(ind,:)-gamma_mean_trace(ind,:)))/norm(v_truth(ind,:));
%   pause(.2)
  
end