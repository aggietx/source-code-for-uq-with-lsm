disp('enkbf_layer')
close all
MM=50;
K_max=Kmax;
% T=400;dt=1/1000;N=T/dt;L1=2*pi;L=36*2;fprintf('L is %d \n',L);
% x=zeros(L,N);sig_ex=0.25;
% y=x;
% x(:,1)=L1*rand(L,1);
% y(:,1)=L1*rand(L,1);
% H1=1;H2=1/2;
fprintf('MM is %d\n',MM);
fprintf('L is %d\n',L);
fprintf('sig ex is %2.2f\n',sig_ex);
uens=zeros(1,MM);
psikens=zeros(Kmax,MM);
% uens(:,1)=u(1)+0.1*randn(MM,1);
omega1=H1/sqrt(2);omega3=H2/sqrt(2);
    du=0.0125;dv1=0.0125;dv2=0.0125;dv3=0.0125;dv4=0.0125;
     dkrep=repmat((dk),1,MM);
     sigmakm=diag(sigmak);
     estu=zeros(1,T);
     estpsi=zeros(Kmax,T);
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end
 x_loc = [x(:,i-1),y(:,i-1)];%+sig_ex*sqrt(dt)*randn(L,2);  
     diff_x1 = x(:,i) - x(:,i-1); diff_x2 = x(:,i) - x(:,i-1) + L1; diff_x3 = x(:,i) - x(:,i-1) - L1;  
    diff_y1 = y(:,i) - y(:,i-1); diff_y2 = y(:,i) - y(:,i-1) + L1; diff_y3 = y(:,i) - y(:,i-1) - L1;  
    diff_xtemp = min(abs(diff_x1), abs(diff_x2)); diff_x_index = min(abs(diff_x3), diff_xtemp);
    diff_ytemp = min(abs(diff_y1), abs(diff_y2)); diff_y_index = min(abs(diff_y3), diff_ytemp);
    diff_x1_index = (diff_x_index == abs(diff_x1)); diff_x2_index = (diff_x_index == abs(diff_x2)); diff_x3_index = (diff_x_index == abs(diff_x3)); 
    diff_y1_index = (diff_y_index == abs(diff_y1)); diff_y2_index = (diff_y_index == abs(diff_y2)); diff_y3_index = (diff_y_index == abs(diff_y3)); 
    diff_x = diff_x1 .* diff_x1_index + diff_x2 .* diff_x2_index + diff_x3 .* diff_x3_index;
    diff_y = diff_y1 .* diff_y1_index + diff_y2 .* diff_y2_index + diff_y3 .* diff_y3_index;
    diff_xy = [diff_x; diff_y];%%%%dz

    
Gy1d = [zeros(L,1),exp(1i * x_loc *  kk1dfull).*repk1d ];
Gu=[ones(L,1),zeros(L,2*Kmax)];
G1=Gu;G2=(Gu+Gy1d);   
fullunknow=zeros(1+2*Kmax,MM);
fullunknow(1,:)=uens;fullunknow(2:2:end,:)=psikens;
fullunknow(3:2:end,:)=conj(psikens);
esthat=fullunknow;
gx1=G1*fullunknow;gx2=G2*fullunknow;
gx=[gx1;gx2];
xbar=mean(fullunknow,2)*ones(1,MM);%% x bar( x here is actually uhat, z is acutally position)
PG=(fullunknow-xbar)*gx'/(MM-1);%%% x_i minu x bar
dz=(diff_xy*ones(1,MM));
daterm=-PG/2/(sig_ex^2)*(gx*dt-dz+sig_ex*sqrt(dt)*randn(2*L,MM));

if max(abs(imag(gx(:))))>10^(-10)
    disp('error complex vel')
end

psikens=psikens-dkrep.*psikens*dt+...
        1i*(kk0*lx.*( (beta./(kk0.^2)/lx^2)*ones(1,MM) -repmat(uens,Kmax,1))).*psikens*dt+...
        (1i*((kk0*lx./(kk0.^2)/lx^2.*hk)*ones(1,MM)).*repmat(uens,Kmax,1)*dt)+...
        1/sqrt(2)*sigmakm*(randn(Kmax,MM)+1i*randn(Kmax,MM))*sqrt(dt)+daterm(2:2:end,:);    
% psistar=real(psikens)-1i*imag(psikens);
uens=uens-du*uens*dt-1i*lx*sum((kk0.*hk)*ones(1,MM).*conj(psikens))*dt+...
        -1i*lx*sum(-(kk0.*hkstar)*ones(1,MM).*psikens)*dt+sigmau*randn(1,MM)*sqrt(dt)+daterm(1,:);
    if mod(i,1/dt)==0
      estu(i*dt)=mean(uens);  
        estpsi(:,i*dt)=mean(psikens,2);
    end
end

rmscc(estu,u(1/dt:1/dt:end)',1);
% rmscc(estpsi,psik(:,1/dt:1/dt:end),1);
for i=1:2;1:10;
rmscc(estpsi(i,:),psik(i,1/dt:1/dt:end),1);
end
estuhat=zeros(1+2*Kmax,T);
estuhat(1,:)=estu;
estuhat(2:2:end,:)=estpsi;
estuhat(3:2:end,:)=conj(estpsi);
save estuhat estuhat;save uhate uhate
uexact=u(1/dt:1/dt:end); psikexact=psik(:,1/dt:1/dt:end);
 compute_errorlayer1d
  errest=[mean(rmstime(10:end)),mean(rmsall),mean(ccall)];
  fprintf('test case is %d\n',testcase);
% clear x y u  psik v1ens v2ens v3ens v4ens u uass psi1  estuhat2  uens psi1ens psi2  psi2ens 