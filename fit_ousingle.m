% clear;
close all
% Set rng seed for reproducability.
rng(1)

% stepSize = 0.01; % Step size for numerical integration.
stepSize = 1/24; % Step size for numerical integration.
tFinal = 365;
tArray = 0:stepSize:tFinal;
nIterations = length(tArray);

% Setting true parameters.
dTruth = 0.5;3; % damping
phiTruth = 0; % imaginary component of damping
fTruth = 0.0901 - 0.0434*1i;2; % forcing
sigTruth =0.2; 2; % stochastic forcing

% Initializing with NaN helps with debugging.
uArray = nan(1, nIterations);
uArray(1) = 0.1061;

% Euler-Maruyama scheme.
for iIteration = 1:nIterations-1
    u = uArray(iIteration);
    
    % I prefer short variable names when implementing formulae.
    d = dTruth;
    phi = phiTruth;
    f = fTruth;
    sig = sigTruth;
    
    dt = stepSize;
    dW = sqrt(dt)*randn(1);
    
    du = ((-d + phi*1i)*u + f)*dt + (sig+1i*sig)*dW;
    
    uArray(iIteration+1) = u + du;
end
 ix=-1; iy=-1;   
[a]=find(kk(1,:)==ix);
[b]=find(kk(2,:)==iy);
ind=intersect(a,b);
% dt=0.01;load u_hat;uArray=u_hat(ind,:);
dt=5.0000e-04;load u_hatsample;uArray=u_hatsample;
funfit=0;
naut=800;%%% ## of points to compute ACF
Nf = naut; xft = 0:dt:naut*dt; % number of lags
    fitfun = fittype('exp(-b*x)*cos(c*x)','coeff', {'b','c'});
    fitfuni = fittype('-exp(-b*x)*sin(c*x)','coeff', {'b','c'});
    fitfunii = fittype('exp(-b*x)*sin(c*x)','coeff', {'b','c'});

sized=size(uArray,2);

streamf=uArray;
Ek=var(streamf,1);
meank=mean(streamf);
% autok=zeros(sized,1);
% for tao=0:sized-1
% autok(tao+1)=sum((streamf(1:end-tao)-meank).*conj(streamf(tao+1:end)-meank));
% end
autok=xcorr(streamf-meank);
autok=autok/Ek/sized;
% autok=autok(sized:-1:1);
autok=autok(sized:end);
% % % sumautok=sum(autok(1:naut))*dt;
sumautok=trapz(autok(1:naut))*dt;
acf=autok;
if funfit
fitr = fit(xft',real(acf(1:Nf+1))',fitfun,'Lower', [0,0], 'Upper', [1e5,pi], 'StartPoint',[1 1]);
Tk= integrate(fitr, 1e5, 0); % acfr = fitr(xft); trapz(xft,acfr);
if(imag(acf(2))<0)
fiti = fit(xft',imag(acf(1:Nf+1))',fitfuni,'Lower', [0,0], 'Upper', [1e5,pi],'StartPoint',[1 1]);
thetak = -integrate(fiti, 1e5, 0);
else
fiti = fit(xft',imag(acf(1:Nf+1))',fitfunii,'Lower', [0,0], 'Upper', [1e5,pi],'StartPoint',[1 1]);
thetak = -integrate(fiti, 1e5, 0);
end
else
Tk=real(sumautok);thetak=-imag(sumautok);
end
dk=Tk/(Tk^2+thetak^2);
omegak=-thetak/(Tk^2+thetak^2);
sigmak2=2*Ek*dk;
% mu(ind)=meank*dk(ind);
mu=meank/(Tk+sqrt(-1)*thetak);
% fprintf('exact are %2.2f %2.2f,%2.2f,%2.3f\n',dTruth,phiTruth,fTruth,sigTruth)
% fprintf('exact are %2.2f %2.2f,%2.2f,%2.3f\n',-L_u(1),0,F_u(1),1 / sqrt(2) * sigma_B)
fprintf('fit are %2.2f %2.2f,%2.2f,%2.3f\n',dk,omegak,mu,sqrt(sigmak2)/sqrt(2))