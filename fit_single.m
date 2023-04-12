clear;
close all
rng(1)
dt=1/5000;
% load u_hatsample;u_hatsample=u_hatsample(1/dt:end);plot(dt:dt:dt*length(u_hatsample),real(u_hatsample));%%%gr
% load u_hatsample1;plot(dt:dt:dt*length(u_hatsample1),real(u_hatsample1))%%%gb
% uArray=u_hatsample;
% uArray=u_hatsample1;
id=4;
load u_hat;uArray=u_hat(id,1:end);
% load cg_m;uArray=cg_m(id,:);
plotpdf(real(uArray));%plotacf1(real(uArray),100/dt,dt);
figure();plot(dt:dt:dt*length(uArray),real(uArray));setgca(18)%%%gb
funfit=0;
naut=floor(length(uArray)/1.2);100000;%%% ## of points to compute ACF
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
figure(); plot(real(autok));setgca(18)
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
% fprintf('fit are %2.2f %2.2f,%2.2f,%2.3f\n',dk,omegak,mu,sqrt(sigmak2)/sqrt(2))
fprintf('fit are %2.2f %2.2f,%2.2f,%2.3f\n',dk,omegak,mu,sqrt(abs(sigmak2))/sqrt(2))