function [sig,dis]=computeRentropy(u_hat,gamma_mean_trace,gamma_cov_trace,Rexact,diagcov)

N1=size(u_hat,2);
dimuhat=size(u_hat,1);

%% filter vs true, exact parameters
sig=zeros(N1,1);dis=sig;

for i=2:N1
mu1=u_hat(:,i);
mu2=gamma_mean_trace(:,i);
if diagcov~=1
Rfil=gamma_cov_trace(:,:,i);
else
    Rfil=diag(gamma_cov_trace(:,i));
end
% invR=inv(Rexact);
invR=Rexact\eye(dimuhat);
sig(i)=(mu1-mu2)'*invR*(mu1-mu2);
dis(i)=trace(Rfil*invR)-dimuhat-log(det(Rfil*invR));
sig(i)=sig(i)/2;dis(i)=dis(i)/2;
if max(abs(imag(sig(i))))>10^(-7) || max(abs(imag(dis(i))))>10^(-7)
    disp('error, complex')
else
    sig(i)=real(sig(i));
    dis(i)=real(dis(i));
end
end
dis(1)=dis(2);
sig(1)=sig(2);
fprintf('mean sig, disp, sum are (filter vs true) %2.2f %2.2f %2.2f\n',mean(sig),mean(dis),mean(sig)+mean(dis));

