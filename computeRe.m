gap=100;N1=N/gap;
Rexact=zeros(dimuhat,dimuhat);
for j=1:dimuhat
streamf=u_hat(j,:);
Rexact(j,j)=var(streamf,1);
end
Rexact(1,1)=1;
sig=zeros(N1,1);dis=sig;
%% true vs filter, exact parameters
for i=gap:gap:N
mu1=u_hat(:,i);
mu2=exact_gamma_mean_trace(:,i);
% R=exact_gamma_cov_trace(:,i);R(1)=1;invR=diag(1./(R));
R=exact_gamma_cov_trace(:,:,i);R(1)=1;invR=inv(R);
sig(i/gap)=(mu1-mu2)'*invR*(mu1-mu2);
dis(i/gap)=trace(Rexact*invR)-dimuhat-log(det(Rexact*invR));
if max(abs(imag(sig(i/gap))))>10^(-7) || max(abs(imag(dis(i/gap))))>10^(-7)
    disp('error, complex')
else
    sig(i/gap)=real(sig(i/gap));
    dis(i/gap)=real(dis(i/gap));
end
end

fprintf('mean sig and disp are (true vs filter) %2.3f %2.3f\n',mean(sig),mean(dis));


%% true vs filter, estimated parameters

sig1=zeros(N1,1);dis1=sig1;
for i=gap:gap:N
mu1=u_hat(:,i);
mu2=gamma_mean_trace(:,i);
% R=exact_gamma_cov_trace(:,i);R(1)=1;invR=diag(1./(R));
R=gamma_cov_trace(:,:,i);R(1)=1;invR=inv(R);
sig1(i/gap)=(mu1-mu2)'*invR*(mu1-mu2);
dis1(i/gap)=trace(Rexact*invR)-dimuhat-log(det(Rexact*invR));
if max(abs(imag(sig1(i/gap))))>10^(-7) || max(abs(imag(dis1(i/gap))))>10^(-7)
    disp('error, complex')
else
    sig1(i/gap)=real(sig1(i/gap));
    dis1(i/gap)=real(dis1(i/gap));
end
end

fprintf('mean  sig and disp are (true vs filter) %2.3f %2.3f\n',mean(sig1),mean(dis1));

%% filter vs true, exact parameters
sig2=zeros(N1,1);dis2=sig2;

for i=gap:gap:N
mu1=u_hat(:,i);
mu2=exact_gamma_mean_trace(:,i);
% R=exact_gamma_cov_trace(:,i);R(1)=1;invR=diag(1./(R));
% R=exact_gamma_cov_trace(:,:,i);R(1)=1;invR=inv(R);
Rfil=exact_gamma_cov_trace(:,:,i);
invR=inv(Rexact);
sig2(i/gap)=(mu1-mu2)'*invR*(mu1-mu2);
dis2(i/gap)=trace(Rfil*invR)-dimuhat-log(det(Rfil*invR));
if max(abs(imag(sig2(i/gap))))>10^(-7) || max(abs(imag(dis2(i/gap))))>10^(-7)
    disp('error, complex')
else
    sig2(i/gap)=real(sig2(i/gap));
    dis2(i/gap)=real(dis2(i/gap));
end
end

fprintf('mean sig and disp are (filter vs true) %2.3f %2.3f\n',mean(sig2),mean(dis2));

%% filter vs true, estimate parameters
sig3=zeros(N1,1);dis3=sig3;

for i=gap:gap:N
mu1=u_hat(:,i);
mu2=exact_gamma_mean_trace(:,i);
% R=exact_gamma_cov_trace(:,i);R(1)=1;invR=diag(1./(R));
% R=exact_gamma_cov_trace(:,:,i);R(1)=1;invR=inv(R);
Rfil=gamma_cov_trace(:,:,i);
invR=inv(Rexact);
sig3(i/gap)=(mu1-mu2)'*invR*(mu1-mu2);
dis3(i/gap)=trace(Rfil*invR)-dimuhat-log(det(Rfil*invR));
if max(abs(imag(sig3(i/gap))))>10^(-7) || max(abs(imag(dis3(i/gap))))>10^(-7)
    disp('error, complex')
else
    sig3(i/gap)=real(sig3(i/gap));
    dis3(i/gap)=real(dis3(i/gap));
end
end

fprintf('mean sig and disp are (filter vs true) %2.3f %2.3f\n',mean(sig3),mean(dis3));


