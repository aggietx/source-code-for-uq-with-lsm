%%% should run layered_cg first 
load Delta_vx_real;load sigma_wx_var_real;load xxvx_real
load Delta_vx_imag;load sigma_wx_var_imag;load xxvx_imag

MM=10;N=1/dt;
psikens=zeros(Kmax,N,MM);
uens=zeros(1,N,MM);
xens=zeros(L,N,MM);
yens=zeros(L,N,MM);
initialpsik=psik0(:,end);
initialu=u0(:,end);
rng(10);
disp('prediction with layered model....')
xinit=randn(L,1)*pi/10+pi;
yinit=randn(L,1)*pi/10+pi;
for j=1:MM
    psik=zeros(Kmax,N);u=zeros(1,N);
temppsik=initialpsik;%% initial condition
tempu=initialu;
  psik(:,1)=initialpsik;u(1)=initialu;
x=zeros(L,N);y=zeros(L,N);
% x(:,1)=x0(:,end);y(:,1)=y0(:,end);
x(:,1)=xinit;y(:,1)=yinit;
fprintf('rest ensemble is %d\n',MM-j);
for i=2:N
%         if mod(i,floor(N/5)) == 0
% %         disp(i*dt)
%  fprintf('rest step is %d\n',N-i);
%         end
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
    
    if max(abs(real(psik(:,i))))>1000
        disp('error,blow up');
        i
        return
    end
end
xens(:,:,j)=x;yens(:,:,j)=y;
psikens(:,:,j)=psik;uens(:,:,j)=u;
% disp('.......')
end


disp('prediction with ou....')
psikens2=zeros(Kmax,N,MM);
uens2=zeros(1,N,MM);
xens2=zeros(L,N,MM);
yens2=zeros(L,N,MM);
dim=length(dkfit);
initialpsik=gamma_mean_trace(2:2:end,end);
initialu=gamma_mean_trace(1,end);
for j=1:MM
    fprintf('rest ensemble is %d\n',MM-j);
psikens2(:,1,j)=initialpsik;uens2(:,1,j)=initialu;

ufree=zeros(dim, N);
uold=zeros(dim,1);
uold(2:2:end)=initialpsik;uold(1)=initialu;
uold(3:2:end)=conj(initialpsik);
x=zeros(L,N);y=zeros(L,N);
% x(:,1)=x0(:,end);y(:,1)=y0(:,end);
x(:,1)=xinit;y(:,1)=yinit;
   for i=2:N 
       ureal=real(uold([1,2:2:end]));uimag=imag(uold([3:2:end]));
     sigma_Z_index = round( (ureal - xxvx_real)./Delta_vx_real');
     sigma_real=zeros(1+Kmax,1);
     for ii=1:1+Kmax
    if sigma_Z_index(ii) <= 0 || sigma_Z_index(ii) > length(sigma_wx_var_real)
        sigma_real(ii) = 0;
    else
        sigma_real(ii) = sigma_wx_var_real(ii,sigma_Z_index(ii));
    end      
     end
  
     sigma_Z_index = round( (uimag - xxvx_imag(2:end))./Delta_vx_imag(2:end)');
     sigma_imag=zeros(1+Kmax,1);
     for ii=2:1+Kmax
    if sigma_Z_index(ii-1) <= 0 || sigma_Z_index(ii-1) > length(sigma_wx_var_imag)
        sigma_imag(ii) = 0;
    else
        sigma_imag(ii) = sigma_wx_var_imag(ii,sigma_Z_index(ii-1));
    end      
     end
  

  idrange=[1,2:2:2*Kmax+1];
  tempu=uold(idrange)+(-dkfit(idrange)+1i*omegakfit(idrange)).*uold(idrange)*dt+...
      Ffit(idrange)*dt+sqrt(dt)*(sigma_real.*randn(Kmax+1,1)+sigma_imag.*randn(Kmax+1,1));
u1=zeros(1+2*Kmax,1);u1(idrange)=tempu;u1(3:2:2*Kmax+1)=conj(tempu(2:end));
  % %    u1=uold+(-dkfit+1i*omegakfit).*uold*dt+Ffit*dt+sqrt(dt)*sigmafit*randn(dim,1);
  
   
    G0=zeros(L,2*Kmax+1);G0(:,1)=1;
    G1 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (-ones(L,1) * transpose(dY));
    G2 = exp(1i * x_loc * kk1dfullf*2*pi/L1).* (ones(L,1) * transpose(dX));
    G1=G1+G0;G2=G2+G0;
    ux=G1*uold;uy=G2*uold;
    

if max(abs(imag(ux(:))))>10^(-6) || max(abs(imag(uy(:))))>10^(-6)
    disp('error, complex velicity')
end
ux=real(ux);uy=real(uy);
 uold=u1;
%% update position
    x(:,i) = x(:,i-1) + ux * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in x
    y(:,i) = y(:,i-1) + uy * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in y
    x(:,i) = mod(x(:,i),L1);
    y(:,i) = mod(y(:,i),L1); 
     ufree(:,i)=u1;

    end
xens2(:,:,j)=x;yens2(:,:,j)=y;
psikens2(:,:,j)=ufree(2:2:end,:);uens2(:,:,j)=ufree(1,:);
% disp('.......')   
end
save psikens psikens
save uens uens
save psikens1 psikens2
save uens1 uens
a=psikens(1,:,:);a=reshape(a(:),size(a,2),size(a,3));
figure();plot(real(a))
a=xens(1,:,:);a=reshape(a(:),size(a,2),size(a,3));
figure();plot(real(a))

save xens xens;save yens yens;save xens2 xens2;save yens2 yens2

% get xens.mat yens.mat xens1.mat yens1.mat
return
load xens;load yens;
load xens2;load yens2;

% close all
figure();
for id=1:16

range=1:800;
subplot(4,4,id)
for iens=1:10
    xtemp=fixtraj(xens(id,range,iens));
    ytemp=fixtraj(yens(id,range,iens));
plot(xtemp,ytemp,'k');hold on
end
plot(xtemp(1),ytemp(1),'or','linewidth',2);hold on
% xlim([0,2*pi]);ylim([0,2*pi])
 xlim([-12,12]);ylim([-12,12])
% figure();

for iens=1:10
    xtemp=fixtraj(xens2(id,range,iens));
    ytemp=fixtraj(yens2(id,range,iens));
plot(xtemp,ytemp,'b');hold on
end
xlim([-12,12]);ylim([-12,12])
setgca(16)
end