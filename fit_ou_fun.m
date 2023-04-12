function [dk,omegak,mu,sigmak2]=fit_ou_fun(funfit,K_max,maxlag,dt,u_hat,kk)
% [dk,omegak,sigmak2,mu]=fit_ou_fun(0,K_max,500,dt,u_hat,kk);
% Kmax=K_max;
% Ktol=18;%funfit=0;
Nf = maxlag; xft = 0:dt:maxlag*dt; % number of lags
    fitfun = fittype('exp(-b*x)*cos(c*x)','coeff', {'b','c'});
    fitfuni = fittype('-exp(-b*x)*sin(c*x)','coeff', {'b','c'});
    fitfunii = fittype('exp(-b*x)*sin(c*x)','coeff', {'b','c'});
    dimuhat=size(u_hat,1);
dk=zeros( dimuhat,1);
omegak=zeros( dimuhat,1);
sigmak2=zeros( dimuhat,1);
mu=zeros( dimuhat,1);
ts=1;
sized=size(u_hat(:,ts:end),2);
tic
% % for ix=-Kmax:Kmax
% %     for iy=-Kmax:Kmax
% for ix=-10:10;-ny/2+1:ny/2;
%     for iy=-10:10;-nx/2+1:nx/2;        
% [a]=find(kk(1,:)==ix);
% [b]=find(kk(2,:)==iy);
% ind=intersect(a,b);
u_hat=u_hat(:,ts:end);
% parfor ind=1:size(dk,1)
for ind=1:size(dk,1)    
streamf=u_hat(ind,:);

Ek=var(streamf,1);
if Ek>10^(-30) %%%&& sqrt(ix^2+iy^2)<Ktol
meank=mean(streamf);
autok=xcorr(streamf-meank);
autok=autok/Ek/sized;
autok=autok(sized:end);
sumautok=trapz(autok(1:maxlag))*dt;
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
dk(ind)=Tk/(Tk^2+thetak^2);
omegak(ind)=-thetak/(Tk^2+thetak^2);
sigmak2(ind)=2*Ek*dk(ind);
% mu(ind)=meank*dk(ind);
mu(ind)=meank/(Tk+sqrt(-1)*thetak);
% mu(ind)=real(meank)*dk(ind)+sqrt(-1)*imag(meank)*omegak(ind);
% ix
% iy
% sigmak2(ind)
%         fname = sprintf('fstreampsi%d_%d_%d.bin', nx,ix,iy);
%      fid=fopen(fname,'w+');fwrite(fid,datap,'single');fclose(fid);%%%
end
%     end
% end

end
sigmak2temp=sigmak2;
sigmak2=sqrt(abs(sigmak2))/sqrt(2);
[a]=find(kk(1,:)==0);
[b]=find(kk(2,:)==0);
ind=intersect(a,b);%% mode (0,0), just take square root
for ii=1:length(ind)
sigmak2(ind(ii))=sqrt(sigmak2temp(ind(ii)));
end