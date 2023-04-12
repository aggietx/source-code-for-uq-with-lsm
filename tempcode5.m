
Nit=1;maxlag=5000;funfit=0;tfit=1;fprintf('lag si %d \n',maxlag);fprintf('total iteration %d \n',Nit);
fprintf('funfit is %d \n',funfit);
estway=1;%%0:filter 1:smoother
sig_ex1=sig_ex*1.00;%%% use inflation;
%% initial guess
K_max=10;Kuse=K_max;if K_max~=Kuse;reduced=1;else reduced=0;end;
fprintf('Kmax and Kuse are %d %d\n',K_max,Kuse);
if Kuse>K_max;Kuse=K_max;reduced=0;end
dimuhat=(2*K_max+1);Kmax=K_max;

% ix=min(K_max,k_0);iy=min(K_max,k_0);ind=getind(ix,iy,kk);%% sample plot

redind=1:dimuhat;
u_hat=zeros(1+2*Kmax,N);
u_hat(1,:)=u;u_hat(2:2:end,:)=psik;
u_hat(3:2:end,:)=conj(psik);

[dke,omegake,mue,sigmak2e]=fit_ou_fun(funfit,Kmax,maxlag,dt,u_hat(:,tfit/dt:end),kk1dfullf);
sigmafit=form_sigmatrix1d(dimuhat,sigmak2e,kk1dfullf,Kmax);
N1=N;1000/dt;
usim=zeros(1+2*Kmax,N1);
for i=2:N1
    usim(:,i)=usim(:,i-1)+(-dke+1i*omegake).*usim(:,i-1)*dt+dt*mue+...
        sqrt(dt)*sigmafit*randn(1+2*Kmax,1);
end

% return
% get a.mat b.mat a1.mat b1.mat

% [a1,b1]=plotpdf(usim(1,:));save a1 a1;save b1 b1
% 
% [a,b]=plotpdf(u);save a a;save b b;

[a1,b1]=plotpdf(real(usim(4,:)));save a1 a1;save b1 b1

[a,b]=plotpdf(real(psik(2,:)));save a a;save b b;

return
close all;
load a;load b;load a1;load b1
figure();
plot(b,a,'r');hold on;
plot(b1,a1,'b');

i=1;
[a1,b1]=plotpdf(usim(2*i,:));
[a,b]=plotpdf(psik(i,:));