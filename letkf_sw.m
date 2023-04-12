%%%% letkf for shallow water
clear;close all;rng(1);
dt=1/200;T=20;N=floor(T/dt);regularizeval=10^(-3);
K_max=1;Kuse=1;L=12*3;diagcov=0;if K_max~=Kuse;reduced=1;else reduced=0;end;
if Kuse>K_max;Kuse=K_max;reduced=0;end
eps=1;E_0=1;E_1=0.12;k_0=2;alpha=3;5/3;nosmoother=1;backsample=0;
fg=0;fr=0;
dimuhat0=(2*K_max+1)^2;dimuhat=dimuhat0*3;Kmax=K_max;
[ky,kx]=meshgrid([0:K_max,-K_max:-1],[0:K_max,-K_max:-1]);
kk=[kx(:),ky(:)]';

rk1 = [1./sqrt(kk(1,:).^2 + kk(2,:).^2+1) .* (-1i * kk(2,:));
    1./sqrt(kk(1,:).^2 + kk(2,:).^2+1) .* (1i * kk(1,:))];
rk1(:,1) = [0;0]; 
rk2 = [1./sqrt(kk(1,:).^2 + kk(2,:).^2)/sqrt(2)./sqrt(kk(1,:).^2 + kk(2,:).^2 + 1) .* (1i * kk(2,:) + kk(1,:) .* sqrt(kk(1,:).^2 + kk(2,:).^2 + 1));
    1./sqrt(kk(1,:).^2 + kk(2,:).^2)/sqrt(2)./sqrt(kk(1,:).^2 + kk(2,:).^2 + 1) .* (-1i * kk(1,:) + kk(2,:) .* sqrt(kk(1,:).^2 + kk(2,:).^2 + 1))];
% rk2(:,1) =[1i;1]/sqrt(2);
rk2(:,1) =[0;0]/sqrt(2);
rk3 = -[1./sqrt(kk(1,:).^2 + kk(2,:).^2)/sqrt(2)./sqrt(kk(1,:).^2 + kk(2,:).^2 + 1) .* (1i * kk(2,:) - kk(1,:) .* sqrt(kk(1,:).^2 + kk(2,:).^2 + 1));
    1./sqrt(kk(1,:).^2 + kk(2,:).^2)/sqrt(2)./sqrt(kk(1,:).^2 + kk(2,:).^2 + 1) .* (-1i * kk(1,:) - kk(2,:) .* sqrt(kk(1,:).^2 + kk(2,:).^2 + 1))];
% rk3(:,1) = [-1i;1]/sqrt(2);
rk3(:,1) = [0;0]/sqrt(2);

redind=getreducedind(kk,K_max,Kuse);redind=[redind;redind+dimuhat0;redind+2*dimuhat0];
kk0=kk;kk=[kk,kk,kk];negindkk=zeros(dimuhat0,1);
for i=1:length(kk0)
ix=kk0(1,i);iy=kk0(2,i);[a]=find(kk0(1,:)==-ix);[b]=find(kk0(2,:)==-iy);negindkk(i)=intersect(a,b); 
end
% ix=min(K_max,k_0);iy=min(K_max,k_0);ind=getind(ix,iy,kk);%% sample plot
ix=K_max;iy=K_max;ind=getind(ix,iy,kk);%% sample plot
rk=zeros(2,dimuhat);rk(:,1:2*dimuhat0)=[rk1,rk2];
rk(:,negindkk+2*dimuhat0)=rk3;%%% guarantee the symmetry
if reduced==0;redind=(1:dimuhat)';end;
rkred=rk(:,redind);kkred=kk(:,redind);
L1=2*pi;%%% size of the domain
xexact=zeros(L,N);
yexact=xexact;

sig_ex=0.0;%%% noise coefficient in tracer equation
u_hat=zeros(dimuhat,N);
% % dkexact=-0.5*ones(dimuhat,1);
% % omegaexact0=0.5;omegaexact=formomega(dimuhat,omegaexact0,kk);
% % Fuexact=0.2*ones(dimuhat,1);
% % sigmauexact=0.5*ones(dimuhat,1);
% % sigmauhat=form_sigmatrix1(dimuhat,sigmauexact,kk);

[dkexactb,omegaexactb,Fuexactb,sigmauhatb,eneb,sigmak2exactb]=form_coeff(dimuhat0,kk0,E_0,k_0,alpha,Kmax,fg);%%%%gb
[dkexactg,omegaexactg,Fuexactg,sigmauhatg,eneg,sigmak2exactg]=form_coeff(dimuhat0,kk0,E_1,k_0,alpha,Kmax,fr);%%%%gravity
dkexact=-[dkexactb;dkexactg;dkexactg];
omegaexact=[zeros(dimuhat0,1);1/eps*sqrt(kk0(1,:).^2 + kk0(2,:).^2)';-1/eps*sqrt(kk0(1,negindkk).^2 + kk0(2,negindkk).^2)'];
Fuexact=[Fuexactb;Fuexactg;Fuexactg];
sigmauhat=zeros(dimuhat,dimuhat);sigmauhat(1:dimuhat0,1:dimuhat0)=sigmauhatb;
sigmauhat(1+dimuhat0:end,1+dimuhat0:end)=form_gravity_sigma(kk0,sigmak2exactg);
% return
% % % sigmak2exact=sigmak2exact/sqrt(2);
% fixedcov=diag((sigmak2exact.^2)./(-dkexact+sqrt(dkexact.^2+L/sig_ex/sig_ex.*sigmak2exact.^2)));
% fixedcov(1,1)=1;fixedcov=fixedcov(redind,redind);
%% generate true signal
% % utest=zeros(size(u_hat));utest1=zeros(size(u_hat
disp('generating true signal...')
for i=2:N
    if mod(i,floor(N/5)) == 0;fprintf('rest step is %d\n',N-i);end
    u_hat(:,i)=u_hat(:,i-1)+(dkexact+1i*omegaexact).*u_hat(:,i-1)*dt+...
        Fuexact*dt+1*sqrt(dt)*sigmauhat*randn(dimuhat,1);
    if max(abs(real(u_hat(:,i))))>10^8
        disp('error, data blow up')
    end
    tempcode=u_hat(:,i);
% %     utest1(:,N-i+1)=tempcode;
% %     utest(:,i)=tempcode;
    
end
xexact(:,1)=L1*rand(L,1);
yexact(:,1)=L1*rand(L,1);
Dim_Grid = 40;
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
uexact=zeros(Dim_Grid^2,T);vexact=zeros(Dim_Grid^2,T);
uexactfilter=zeros(Dim_Grid^2,T);vexactfilter=zeros(Dim_Grid^2,T);
uexactsmoother=zeros(Dim_Grid^2,T);vexactsmoother=zeros(Dim_Grid^2,T);
uexactbackward=zeros(Dim_Grid^2,T);vexactbackward=zeros(Dim_Grid^2,T);
    Ga = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(1,:))); % Fourier bases for u
    Gb = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(2,:))); % Fourier bases for v
% return
disp('generating obsveration...')
time=1;
for i=2:N
    x_loc = [xexact(:,i-1),yexact(:,i-1)];
    G1 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
    G2 = (exp(1i * x_loc * kk*2*pi/L1) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
    u = (G1*  u_hat(:,i-1)); 
    v = (G2*  u_hat(:,i-1)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity'  )
      max(abs(imag(u)))
    end
    u=real(u);
    v=real(v);
    xexact(:,i) = xexact(:,i-1) + u * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in x
    yexact(:,i) = yexact(:,i-1) + v * dt + sqrt(dt) * sig_ex * randn(L,1); % floe equation in y
%     xexact=real(xexact);yexact=real(yexact);
if imag(max(abs(xexact(:,i))))>0 || imag(max(abs(yexact(:,i))))>0
    disp('error, complex')
end
    xexact(:,i) = mod(xexact(:,i),L1);
    yexact(:,i) = mod(yexact(:,i),L1);
    if mod(i,1/dt)==0
    u = (Ga*  u_hat(:,i-1)); 
    v = (Gb*  u_hat(:,i-1)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity, record exact'  )
      max(abs(imag(u)))
    end
    uexact(:,time)=real(u);
    vexact(:,time)=real(v);
     
        time=time+1;  
    end
end




