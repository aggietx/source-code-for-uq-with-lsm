%% analytical ene
exactene=ene;
F1=real(allFfit(:,Nit));F2=imag(allFfit(:,Nit));%%%% use estimated parameter to compute the energy
estene=1/2*((F2.^2+ F1.^2)./alldkfit(:,Nit) +allsigmak2(:,Nit).^2*2./(2*alldkfit(:,Nit)) );

% F1=real(allFfit(:,end));F2=imag(allFfit(:,end));%%%% use estimated parameter to compute the energy
% exene=1/2*(allsigmak2(:,end-1).^2*2./(-2*alldkfit(:,end-1)) );
% [a]=computeene1D(Kmax,estene,kk);
% [b]=computeene1D(Kmax,exene,kk);
[exactenekmean,exacteneksum]=computeene1D(Kmax,exactene,kk);
[estenekmean,esteneksum]=computeene1D(Kmax,estene,kk);
% plot2line(exactenekmean,estenekmean);
%% discrete ene
a=real(u_hat);b=imag(u_hat);
exactene1=sum(a.^2+b.^2,2);
% a=real( exact_gamma_mean_trace_smoother);b=imag( exact_gamma_mean_trace_smoother);
a=real( exact_gamma_mean_trace);b=imag( exact_gamma_mean_trace);
smene1=sum(a.^2+b.^2,2);
% a=real( gamma_mean_trace_smoother);b=imag(gamma_mean_trace_smoother);
a=real( gamma_mean_trace);b=imag(gamma_mean_trace);
estene1=sum(a.^2+b.^2,2);

[exactenekmean1,exacteneksum1]=computeene1D(Kmax,exactene1,kk);
[smenekmean1,smeneksum1]=computeene1D(Kmax,smene1,kk);
[estenekmean1,esteneksum1]=computeene1D(Kmax,estene1,kk);
% % % plot2line(exactenekmean1,smenekmean1);
% % % plot2line(smenekmean1,estenekmean1);

sampletraj=zeros(dimuhat,N);
 for i=1:N
    if size(gamma_cov_trace,3)==1
   
    tempR=gamma_cov_trace(:,i);
% sampletraj(:,i)=gamma_mean_trace_smoother(:,i)+sqrt(tempR).*randn(dimuhat,1);
 sampletraj(:,i)=gamma_mean_trace(:,i)+sqrt(tempR).*randn(dimuhat,1);   
     else
       tempR=gamma_cov_trace(:,:,i);
       tempR=diag(tempR);
% sampletraj(:,i)=gamma_mean_trace_smoother(:,i)+sqrt(tempR).*randn(dimuhat,1); 
sampletraj(:,i)=gamma_mean_trace(:,i)+sqrt(tempR).*randn(dimuhat,1); 
    end
end
a=real(sampletraj);b=imag(sampletraj);
sampletrajene1=sum(a.^2+b.^2,2);
[sampleenekmean1,sampleeneksum1]=computeene1D(Kmax,sampletrajene1,kk);
sampletraj=sampletraj(:,gap:gap:end);

% figure();plot(1:length(exactenekmean1),exactenekmean1,'-*r','linewidth',2);hold on;
% plot(1:length(exactenekmean1),smenekmean1,'-*b','linewidth',2);hold on;
% plot(1:length(exactenekmean1),estenekmean1,'-*g','linewidth',2);
% plot(1:length(exactenekmean1),sampleenekmean1,'-*k','linewidth',2);
% lgnd=legend('Exact','Smoother with exact parameters','Smoother with estimated parameters','Reconstruct');
% set(lgnd,'FontSize',14);
% setgca(16);
%% transform estimation to 2d
dkfit2d=plotest2d(Kmax,kk,alldkfit(:,Nit));
omegafit2d=plotest2d(Kmax,kk,allomegakfit(:,Nit));
Ffit2d=plotest2d(Kmax,kk,allFfit(:,Nit));
sigmafit2d=plotest2d(Kmax,kk,allsigmak2(:,Nit)*sqrt(2));

dkexact2d=plotest2d(Kmax,kk,-dkexact);
% omegaexact2d=plotest2d(Kmax,kk,omegaexact);
sigmaexact2d=plotest2d(Kmax,kk,sigmak2exact);

estpara2d=zeros(2*Kmax+1,2*Kmax+1,4);
estpara2d(:,:,1)=dkfit2d;estpara2d(:,:,2)=omegafit2d;
estpara2d(:,:,3)=Ffit2d;estpara2d(:,:,4)=sigmafit2d;

exactpara2d=zeros(2*Kmax+1,2*Kmax+1,4);
exactpara2d(:,:,1)=dkexact2d;exactpara2d(:,:,2)=0;
exactpara2d(:,:,3)=0;exactpara2d(:,:,4)=sigmaexact2d;

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
% savedtraj=[u_hat(:,gap1:gap1:end);
%     exact_gamma_mean_trace_smoother(:,1:gap1:end);
%     gamma_mean_trace_smoother(:,1:gap1:end)];
savedtraj1=[u_hat(:,gap1:gap1:end);
    exact_gamma_mean_trace(:,gap1:gap1:end);
    gamma_mean_trace(:,gap1:gap1:end)];

%%
gap2=gap;
ue=zeros(dimuhat,N);
ufree=zeros(dimuhat,N);
sigmafit=form_sigmatrix1(size(u_hat,1),allsigmak2(:,Nit),kk,Kmax);
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
disp('free run fixed noise GB error');

[rmstimefree,rmsallfree,ccallfree ]=compute_velerrorGBsw(ue(:,gap2:gap2:end),  ufree(:,gap2:gap2:end),rk,kk,L1,Dim_Grid,3);
errest=[errest,[mean(rmstimefree(3:end));mean(rmsallfree);mean(ccallfree )]];
frdata=[ue(:,gap2:gap2:end);ufree(:,gap2:gap2:end)];