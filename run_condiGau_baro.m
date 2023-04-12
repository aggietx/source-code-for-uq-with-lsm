%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script runs full conditional Gaussian equations for the
% layered topographic barotropic model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontSize', 20)

% set model parameters
kmax = 2; k0 = 1;   % ! here to change the total number of modes
kvec = k0*[1:kmax]';
beta = 1;
H1 = 1*1; H2 = 1*1/2;
theta0 = 0;
lm = 1;
lvec = lm*[cos(theta0),sin(theta0)];

% damping and stochastic forcing parameters
hk = zeros(kmax,1);
hk(1) = H1*(1-1i)/2; hk(2) = H2*(1-1i)/2;
for kk=3:kmax
    hk(kk) = 1/kk^2 *exp(-1i*pi/4)/sqrt(2);
end
d_k = zeros(kmax,1); d_k(:) = 1*.0125; 
d_U = 1*.0125;
% sig_k = zeros(kmax,1); sig_k(:) = 10*1/20/sqrt(2); 
% sig_U = 10*1/20/sqrt(2); 
sig_k = zeros(kmax,1); sig_k(:) = 10*1/20/sqrt(2);
sig_U = 10*1/20/sqrt(2);
f_k = zeros(kmax,1); f_U = 0;

% set tracer parameters
d_T = 0.1; kappa_T = .001;  % dissipations
alpha = 1;                  % mean gradient
% dispersion relations for passive tracer
gamma_T = d_T + kappa_T*kvec.^2;

params = struct('beta',beta,'hk',hk, 'kmax',kmax,'kvec',kvec, ...
    'd_k',d_k,'d_U',d_U,'sig_k',sig_k,'sig_U',sig_U,'gamma_T',gamma_T, ...
    'f_k',f_k,'f_U',f_U, 'lvec',lvec,'lm',lm,'alpha',alpha);
            
% set up model integration params.
dt = 1E-3; Tfin = 500;0; tstep = 10; %500, 100
% -------------------------------------------------------------------------
% initialization for topo. mode model
U0 = .1*randn(1);
omek0 = .1*(1*randn(kmax,1)+1i*randn(kmax,1))/sqrt(2);
Tk0 = .1*(1*randn(kmax,1)+1i*randn(kmax,1))/sqrt(2);
[tout,Uout,dUout, vout,Tout, umout,Ruout,dfout,dgout, energy,enstrophy] = ...
          BaroLayered_tracer_condiGau(U0,omek0,Tk0, dt,Tfin,tstep,params);
   
      
um0 = umout(:,end);
Ru0 = Ruout(:,:,end);
% % run only the full dynamical model
% [tout,Uout,dUout, vout,Tout] = BaroLayered_tracer(U0,omek0,Tk0, dt,Tfin,tstep,params);


% -------------------------------------------------------------------------
% plot results
figure
subplot(5,1,1)
plot(tout(1:end),Uout(1,1:end)); ylabel('U');
title('time series of the topographic barotropic flow & tracer solution')
hold on
% subplot(2,1,2)
% plot(tout(1:end),dUout(1,1:end)); ylabel('dU/dt'); xlabel('time')
% figure
subplot(5,1,2)
plot(tout,real(vout(1,:)), '-'); hold on
plot(tout,imag(vout(1,:)), '--');
legend('real','imaginary'); ylabel('v_1');
subplot(5,1,3)
plot(tout,real(vout(2,:)), '-'); hold on
plot(tout,imag(vout(2,:)), '--');
legend('real','imaginary'); ylabel('v_2');
% xlabel('time')

% figure
subplot(5,1,4)
plot(tout,real(Tout(1,:)), '-'); hold on
plot(tout,imag(Tout(1,:)), '--');
legend('real','imaginary'); ylabel('T_1');
title('time series of the passive tracer solution')
subplot(5,1,5)
plot(tout,real(Tout(2,:)), '-'); hold on
plot(tout,imag(Tout(2,:)), '--');
legend('real','imaginary'); ylabel('T_2');
xlabel('time')

figure
subplot(4,1,1)
plot(tout,real(umout(1,:)), '-'); hold on
plot(tout,imag(umout(1,:)), '--');
legend('real','imaginary'); ylabel('m_{v1}');
title('time series of the conditional mean of flow & tracer solution')
subplot(4,1,2)
plot(tout,real(umout(2,:)), '-'); hold on
plot(tout,imag(umout(2,:)), '--');
legend('real','imaginary'); ylabel('m_{v2}');
% xlabel('time')
% figure
subplot(4,1,3)
plot(tout,real(umout(2*kmax+1,:)), '-'); hold on
plot(tout,imag(umout(2*kmax+1,:)), '--');
legend('real','imaginary'); ylabel('m_{T1}');
% title('time series of the conditional mean of tracer solution')
subplot(4,1,4)
plot(tout,real(umout(2*kmax+2,:)), '-'); hold on
plot(tout,imag(umout(2*kmax+2,:)), '--');
legend('real','imaginary'); ylabel('m_{T2}');
xlabel('time')

figure
subplot(4,1,1)
plot(tout,squeeze(Ruout(1,1,:)+Ruout(2*kmax,2*kmax,:)), '-'); hold on
ylabel('r_{v1}');
title('time series of the conditional variance of flow & tracer solution')
subplot(4,1,2)
plot(tout,squeeze(Ruout(2,2,:)+Ruout(2*kmax-1,2*kmax-1,:)), '-'); hold on
ylabel('r_{v2}');
% xlabel('time')
% figure
subplot(4,1,3)
plot(tout,squeeze(Ruout(2*kmax+1,2*kmax+1,:)+Ruout(4*kmax,4*kmax,:)), '-'); hold on
ylabel('r_{T1}');
% title('time series of the conditional variance of tracer solution')
subplot(4,1,4)
plot(tout,squeeze(Ruout(2*kmax+2,2*kmax+2,:)+Ruout(4*kmax-1,4*kmax-1,:)), '-'); hold on
ylabel('r_{T2}');
xlabel('time')