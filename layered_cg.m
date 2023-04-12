%%%% layered model, conditional gaussian
%%%% with direction (lx,ly)
%%%% almost the same fit_layered 
clear;close all;
rng(10);
T=100;dt=1/1000;N=T/dt;L1=2*pi;rk=1;
L=30;Luse=L;min(L,4);Ltotal=L;fprintf('L is %d \n',Luse);
testcase=1;%%%1:non gauss 2:strong non gauss
if testcase==1
    load dkfitc1;load omegakfitc1;load Ffitc1;load sigmafitc1
else
  load dkfitc2;load omegakfitc2;load Ffitc2;load sigmafitc2
end
b1b1t=sigmafit * sigmafit';
fprintf('test case %d\n',testcase);
if testcase==1
constant_cov=0;diagcov=0;
elseif testcase==2
constant_cov=0;diagcov=1;
elseif testcase==3
constant_cov=1;diagcov=1;
end
Kmax=6;p=2;%theta=pi/3;
fprintf('p is %2.2f\n',p);
fprintf('test case is %d\n',testcase);
fprintf('L is %d \n',Luse);
kk0=(1:Kmax)';lx=0.8;ly=sqrt(1-lx^2);fprintf('(lx, ly) are (%2.2f,%2.2f)\n',lx,ly);
dk=0.0125*ones(Kmax,1);du=0.0125;H1=1;H2=1/2;
beta=2;sigmak=zeros(Kmax,1);hk=zeros(Kmax,1);
x=zeros(L,N);sig_ex=.25;
y=x;gap=1;dt0=1*dt;N0=T/dt0;
fprintf('Kmax is %d\n',Kmax);
fprintf('sig ex is %2.2f\n',sig_ex);
x(:,1)=L1*rand(L,1);
y(:,1)=L1*rand(L,1);
if testcase==1;
    sigmau=1/2/sqrt(2);sigmak(1)=1/1/sqrt(2);sigmak(2)=1/2/sqrt(2);
else
    sigmau=1/1/sqrt(2);sigmak(1)=1/2/sqrt(2)/2;sigmak(2)=1/2/sqrt(2)/2;
end
hk(1)=H1/2-H1/2*1i;hk(2)=H2/2-H2/2*1i;

% a=real(psik);b=imag(psik);c=sum(a.^2+b.^2,2);c(1)/c(2)
for ik=3:Kmax
    if testcase==1
    sigmak(ik)=1/ik^p/sqrt(2);
    else
     sigmak(ik)=1/ik^p/sqrt(2)/2;   
    end
    theta=rand(1)*2*pi;
    hk(ik)=sin(theta)/4/(ik)^p-cos(theta)/4/(ik)^p*1i;
    
end
hkstar=real(hk)-1i*imag(hk);
psik=zeros(Kmax,N);u=zeros(1,N);
tempu=u(:,1);
temppsik=psik(:,1);

kk1dfull=zeros(2,2*Kmax);
for i=1:Kmax
   kk1dfull(1,2*i-1)=i*lx;
   kk1dfull(1,2*i)=-i*lx;
   kk1dfull(2,2*i-1)=i*ly;
   kk1dfull(2,2*i)=-i*ly;   
end

kk1dfullf=[[0;0],kk1dfull];
repk1d=repmat(kk1dfull(1,:),L,1)*1i;
kx=kk1dfullf(1,:);ky=kk1dfullf(2,:);
dX = 1i*kx*2*pi/L1;
dY = 1i*ky*2*pi/L1;
dX=dX(:);dY=dY(:);
allvel=zeros(L,N);
rng(10);
% sigmak=zeros(size(sigmak));sigmau=0;
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end
x_loc = [x(:,i-1),y(:,i-1)];          
        %% step1
    k1a=-dk.*temppsik+1i*kk0*lx.*(beta./(kk0.^2)-tempu).*temppsik+...
        1i*kk0*lx./(kk0.^2).*hk.*tempu;
    k1b=-du*tempu-1i*lx*sum(kk0.*hk.* conj(temppsik))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik);
    %% step2
    temppsik2=temppsik+k1a*dt/2; tempu2=tempu+k1b*dt/2;
    k2a=-dk.*temppsik2+1i*kk0*lx.*(beta./(kk0.^2)-tempu2).*temppsik2+...
        1i*kk0*lx./(kk0.^2).*hk.*tempu2;
    k2b=-du*tempu2-1i*lx*sum(kk0.*hk.* conj(temppsik2))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik2);    
    %% step3
    temppsik3=temppsik2+k2a*dt/2; tempu3=tempu2+k2b*dt/2;
    k3a=-dk.*temppsik3+1i*kk0*lx.*(beta./(kk0.^2)-tempu3).*temppsik3+...
        1i*kk0*lx./(kk0.^2).*hk.*tempu3;
    k3b=-du*tempu3-1i*lx*sum(kk0.*hk.*conj(temppsik3))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik3); 
    %% step4
     temppsik4=temppsik3+k3a*dt; tempu4=tempu3+k3b*dt;
     k4a=-dk.*temppsik4+1i*kk0*lx.*(beta./(kk0.^2)-tempu4).*temppsik4+...
        1i*kk0*lx./(kk0.^2).*hk.*tempu4;
     k4b=-du*tempu4-1i*lx*sum(kk0.*hk.* conj(temppsik4))+...
        -1i*lx*sum(-kk0.*hkstar.*temppsik4); 
    ksa=1/6*k1a+1/3*k2a+1/3*k3a+1/6*k4a;
    ksb=1/6*k1b+1/3*k2b+1/3*k3b+1/6*k4b;
psik(:,i)=temppsik+ksa*dt+1/sqrt(2)*sigmak.*(randn(Kmax,1)+1i*randn(Kmax,1))*sqrt(dt); 
% psik(:,i)=temppsik+ksa*dt+(1+1i)/sqrt(2)*sigmak.*randn(Kmax,1)*sqrt(dt); 
% % psik(:,i)=temppsik+ksa*dt+sigmak.*(randn(Kmax,1))*sqrt(dt); 
u(:,i)=tempu+ksb*dt+sigmau*randn*sqrt(dt);
temppsik=psik(:,i);
tempu=u(:,i);
%% get velocity
  
fullunknow=zeros(1+2*Kmax,1);
% fullunknow(1)=u(i-1);
fullunknow(2:2:end)=psik(:,i-1);
fullunknow(3:2:end)=conj(psik(:,i-1));

    sx=-dY.*fullunknow;sy=dX.*fullunknow;
    G = exp(1i * x_loc * kk1dfullf*2*pi/L1);
    ux=G*sx(:)+u(i-1);uy=G*sy(:)+u(i-1);%%% do not forget to add background part
    
%     G0=zeros(L,2*Kmax+1);G0(:,1)=1;
%     G1 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (-ones(L,1) * transpose(dY));
%     G2 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (ones(L,1) * transpose(dX));
%     G1=G1+G0;G2=G2+G0;
%     fullunknow1=fullunknow;fullunknow1(1)=u(i-1);
%     ux1=G1*fullunknow1;uy1=G2*fullunknow1;
%     rnorm(ux1,ux);
%     rnorm(uy1,uy);

if max(abs(imag(ux(:))))>10^(-6) || max(abs(imag(uy(:))))>10^(-6)
    disp('error, complex velicity')
end
ux=real(ux);uy=real(uy);
%% update position
    x(:,i) = x(:,i-1) + ux * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in x
    y(:,i) = y(:,i-1) + uy * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in y
    x(:,i) = mod(x(:,i),L1);
    y(:,i) = mod(y(:,i),L1);    

    allvel(:,i)=sqrt(ux.^2+uy.^2);
    
    if max(abs(real(psik(:,i))))>1000
        disp('error,blow up');
        i
        return
    end
end

%% conditional gaussian
dimuhat=2*Kmax+1;redind=1:dimuhat;nosmoother=1;
Dim_Y=dimuhat;Dim_Yr=length(redind);
Dim_X=2*L;
gamma_mean0 = zeros(Dim_Yr,1);
gamma_cov0 = eye(Dim_Yr)*0.01; 
if constant_cov
    gamma_cov0=fixedcov;
end
% gamma_cov0 =fixedcov;
gamma_mean_trace=zeros(Dim_Y,N);%%% store all
gamma_mean_trace(redind,1) = gamma_mean0;
if nosmoother~=1%%%% only store reduced dof
    if diagcov
gamma_cov_trace=zeros(Dim_Yr,N);
gamma_cov_trace(:,1) = diag(gamma_cov0);        
    else
gamma_cov_trace=zeros(Dim_Yr,Dim_Yr,N);
gamma_cov_trace(:,:,1) = gamma_cov0;
    end
end

A0=zeros(Dim_X,1);
a0=Ffit(redind);
a1=diag(-dkfit(redind)+1i*omegakfit(redind));


invBoB = 1 /sig_ex/ sig_ex * eye(2*L); % inverse of the square of the observational noise
invBoBdiag=1 /sig_ex/ sig_ex * ones(2*L,1);
 G0=zeros(L,2*Kmax+1);G0(:,1)=1;
time=1;
disp('data assimilation (filter)......')
for i = 2:N
    if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
    end
    % observational operator 
    x_loc = [x(:,i-1),y(:,i-1)];
   
    G1 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (-ones(L,1) * transpose(dY));
    G2 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (ones(L,1) * transpose(dX));
    G1=G1+G0;G2=G2+G0;
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
    
% idx=[1:20,1+L:20+L];

    % run the data assimilation for posterior mean and posterior covariance
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * (invBoBdiag.*  (diff_xy - A0*dt-A1 * gamma_mean0 * dt));
if constant_cov
    gamma_cov=fixedcov;
else
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1b1t - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;         
end

    

    % save the posterior statistics
    gamma_mean_trace(redind,i) = gamma_mean;

    % update
    gamma_mean0 = gamma_mean;     
    if diagcov
    gamma_cov0 = diag(diag(gamma_cov));
    else
    gamma_cov0 = gamma_cov;
    end                


% %          if diagcov
% %   gamma_cov_trace(:,i) = diag(gamma_cov0);             
% %          else
% %     gamma_cov_trace(:,:,i) = gamma_cov0;
% %          end
end
rmscc(real( gamma_mean_trace(1,:)),real(u),1);
for i=1:Kmax
rmscc(real( gamma_mean_trace(2*i,:)),real(psik(i,:)),1);
end

N0=N;x0=x;y0=y;
psik0=psik;u0=u;

nx=30;ny=nx;Dim_Grid = nx;L1=2*pi;
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)]; G00=zeros(nx*ny,2*Kmax+1);G00(:,1)=1;
 G1fixed = (exp(1i * xy * kk1dfullf*2*pi/L1) .* (-ones(nx*ny,1) * transpose(dY)))+G00; 
 G2fixed = (exp(1i * xy * kk1dfullf*2*pi/L1) .* (ones(nx*ny,1)* transpose(dX)))+G00;
 
 
uh=zeros(1+2*Kmax,1);
uh(1)=u(end);
uh(2:2:end)=psik(:,end);
uh(3:2:end)=conj(psik(:,end));
 ux1=G1fixed*uh;uy1=G2fixed*uh;
if max(abs(imag( ux1(:))))>10^(-6) || max(abs(imag(uy1(:))))>10^(-6)
    disp('error, complex velicity')
else
    ux1=real(ux1);uy1=real(uy1);
end

 ux2=G1fixed*gamma_mean_trace(:,end);uy2=G2fixed*gamma_mean_trace(:,end);
if max(abs(imag( ux2(:))))>10^(-6) || max(abs(imag(uy2(:))))>10^(-6)
    disp('error, complex velicity')
else
    ux2=real(ux2);uy2=real(uy2);
end
figure();
subplot(2,1,1);
quiver(xx,yy,reshape(ux1,ny,nx),reshape(uy1,ny,nx));
hold on;scatter(x(:,end),y(:,end))
title('exact velocity','fontsize',16);setgca(16)
xlim([0,L1]);ylim([0,L1]);

subplot(2,1,2);
quiver(xx,yy,reshape(ux2,ny,nx),reshape(uy2,ny,nx));
hold on;scatter(x(:,end),y(:,end))
title('filtered velocity','fontsize',16);setgca(16)
xlim([0,L1]);ylim([0,L1]);
%  filenameStrPdf = sprintf('oufitcase%dmode%d.pdf',testcase,i-1);   
print(gcf, '-dpdf', '-r600', 'final_vel.pdf')

return


for i=100:100:N
uh=zeros(1+2*Kmax,1);
uh(1)=u(i);
uh(2:2:end)=psik(:,i);
uh(3:2:end)=conj(psik(:,i));
 ux1=G1fixed*uh;uy1=G2fixed*uh;
if max(abs(imag( ux1(:))))>10^(-6) || max(abs(imag(uy1(:))))>10^(-6)
    disp('error, complex velicity')
else
    ux1=real(ux1);uy1=real(uy1);
end

 ux2=G1fixed*gamma_mean_trace(:,end);uy2=G2fixed*gamma_mean_trace(:,end);
if max(abs(imag( ux2(:))))>10^(-6) || max(abs(imag(uy2(:))))>10^(-6)
    disp('error, complex velicity')
else
    ux2=real(ux2);uy2=real(uy2);
end
% figure();
% subplot(2,1,1);
quiver(xx,yy,reshape(ux1,ny,nx),reshape(uy1,ny,nx));
% hold on;scatter(x(:,end),y(:,end))
% title('exact velocity','fontsize',16);setgca(16)
pause(1)
end