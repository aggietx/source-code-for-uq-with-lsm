
a=real(u_hat);b=imag(u_hat);
exactene1=sum(a.^2+b.^2,2);
egb=exactene1(1:dimuhat0);egb=egb/egb(2);eneb1=eneb/eneb(2);%rnorm(egb,eneb1);
egr=exactene1(1+dimuhat0:2*dimuhat0);egr=egr/egr(2);eneg1=eneg/eneg(2);%rnorm(egr,eneg1);
% sum(eneb)/2/(sum(eneg))
% % sum(exactene1(1:dimuhat0))/2/sum(exactene1(1+dimuhat0:2*dimuhat0))

a=real( exact_gamma_mean_trace_smoother);b=imag( exact_gamma_mean_trace_smoother);
smene1=sum(a.^2+b.^2,2);smene1=smene1/smene1(2);
a=real( gamma_mean_trace_smoother);b=imag(gamma_mean_trace_smoother);
estene1=sum(a.^2+b.^2,2);estene1=estene1/estene1(2);
%% transform estimation to 2d
range=1:dimuhat0;
dkgbfit2d=plotest2d(Kmax,kk0,alldkfit(range,Nit));
omegagbfit2d=plotest2d(Kmax,kk0,allomegakfit(range,Nit));
Fgbfit2d=plotest2d(Kmax,kk0,allFfit(range,Nit));
sigmagbfit2d=plotest2d(Kmax,kk0,allsigmak2(range,Nit)*sqrt(2));

dkgbexact2d=plotest2d(Kmax,kk0,-dkexact(range));
omegagbexact2d=plotest2d(Kmax,kk0,omegaexact(range));
sigmagbexact2d=plotest2d(Kmax,kk0,sigmak2exact(range));

range=1+dimuhat0:2*dimuhat0;
dkgrfit2d=plotest2d(Kmax,kk0,alldkfit(range,Nit));
omegagrfit2d=plotest2d(Kmax,kk0,allomegakfit(range,Nit));
Fgrfit2d=plotest2d(Kmax,kk0,allFfit(range,Nit));
sigmagrfit2d=plotest2d(Kmax,kk0,allsigmak2(range,Nit)*sqrt(2));

dkgrexact2d=plotest2d(Kmax,kk0,-dkexact(range));
omegagrexact2d=plotest2d(Kmax,kk0,omegaexact(range));
sigmagrexact2d=plotest2d(Kmax,kk0,sigmak2exact(range));

estpara2d=zeros(2*Kmax+1,2*Kmax+1,8);
estpara2d(:,:,1)=dkgbfit2d;estpara2d(:,:,2)=omegagbfit2d;
estpara2d(:,:,3)=Fgbfit2d;estpara2d(:,:,4)=sigmagbfit2d;
estpara2d(:,:,5)=dkgrfit2d;estpara2d(:,:,6)=omegagrfit2d;
estpara2d(:,:,7)=Fgrfit2d;estpara2d(:,:,8)=sigmagrfit2d;
exactpara2d=zeros(2*Kmax+1,2*Kmax+1,8);
exactpara2d(:,:,1)=dkgbexact2d;exactpara2d(:,:,2)=omegagbexact2d;
exactpara2d(:,:,3)=0;exactpara2d(:,:,4)=sigmagbexact2d;
exactpara2d(:,:,5)=dkgrexact2d;exactpara2d(:,:,6)=omegagrexact2d;
exactpara2d(:,:,7)=0;exactpara2d(:,:,8)=sigmagrexact2d;
%%% compute the free run using exact and estimated parameter
%%
E1=zeros(N/gap1,1);
E2=zeros(N/gap1,1);
E3=zeros(N/gap1,1);
for ii=1:N/gap1
    i=ii*gap1;
    tempR=Rexact;detR=(det(tempR));%%%% already real
    E1(ii)=dimuhat/2*(1+log(2*pi))+1/2*log(detR);
    
    if size(exact_gamma_cov_trace,3)>1
    tempR=exact_gamma_cov_trace(:,:,i); 
    else
    tempR=diag(exact_gamma_cov_trace(:,i));         
    end
    detR=(det(tempR));
%     if real(detR)/imag(detR)<10^20
%         disp('complex determinant')
%     end
    detR=real(detR);
    E2(ii)=dimuhat/2*(1+log(2*pi))+1/2*log(detR);
        if size(gamma_cov_trace,3)>1
    tempR=gamma_cov_trace(:,:,i); 
        else
    tempR=diag(gamma_cov_trace(:,i));             
        end
       detR=(det(tempR));
%     if real(detR)/imag(detR)<10^20
%         disp('complex determinant')
%     end
    detR=real(detR);
    E3(ii)=dimuhat/2*(1+log(2*pi))+1/2*log(detR);
end
Ent=[E1,E2,E3];

%%
savedtraj=[u_hat(:,gap1:gap1:end);
    exact_gamma_mean_trace_smoother(:,1:gap1:end);
    gamma_mean_trace_smoother(:,1:gap1:end)];

%%
gap2=gap;
ue=zeros(dimuhat,N);
ufree=zeros(dimuhat,N);
% % sigmafit=form_sigmatrix1(size(u_hat,1),allsigmak2(:,Nit),kk,Kmax);
tempsig=allsigmak2(:,Nit);
sigmafit=zeros(dimuhat,dimuhat);
sigmafit(1:dimuhat0,1:dimuhat0)=form_sigmatrix1(dimuhat0,tempsig(1:dimuhat0),kk0,Kmax);
sigmafit(1+dimuhat0:end,1+dimuhat0:end)=form_gravity_sigma(kk0,sqrt(2)*tempsig(1+dimuhat0:2*dimuhat0));

disp('free run with fixed noise (exact parameter and estimated parameter) ..')
for i=2:N
    noise=randn(dimuhat,1);
    ue(:,i)=ue(:,i-1)+(dkexact+1i*omegaexact).*ue(:,i-1)*dt+...
        Fuexact*dt+sqrt(dt)*sigmauhat*noise;
    ufree(:,i)=ufree(:,i-1)+(-alldkfit(:,Nit)+1i*allomegakfit(:,Nit)).*ufree(:,i-1)*dt+...
        allFfit(:,Nit)*dt+sqrt(dt)*sigmafit*noise;
    if max(abs(real(ue(:,i))))>10^8 || max(abs(real(ufree(:,i))))>10^8
        disp('error, data blow up')
    end
%     tempcode=u_hat(:,i);
% %     utest1(:,N-i+1)=tempcode;
% %     utest(:,i)=tempcode;
    
end
disp('free run fixed noise sw error');

 [rmstimefree,rmsallfree,ccallfree ]=compute_velerrorGBsw(ue(:,gap2:gap2:end),  ufree(:,gap2:gap2:end),rk,kk,L1,Dim_Grid,3);
frdata=[ue(:,gap2:gap2:end);ufree(:,gap2:gap2:end)];
errest=[errest,[mean(rmstimefree(3:end));mean(rmsallfree);mean(ccallfree )]];
