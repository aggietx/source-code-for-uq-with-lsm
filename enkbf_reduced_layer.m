disp('enkbf_reduced_layer')
close all
MM=50;K_max=Kmax;
% T=400;dt=1/1000;N=T/dt;L1=2*pi;L=36*2;fprintf('L is %d \n',L);
% x=zeros(L,N);sig_ex=0.25;
% y=x;
% x(:,1)=L1*rand(L,1);
% y(:,1)=L1*rand(L,1);
% H1=1;H2=1/2;
fprintf('MM is %d\n',MM);
fprintf('L is %d\n',L);
fprintf('sig ex is %2.2f\n',sig_ex);
uens=zeros(MM,N);v1ens=uens;v2ens=uens;v3ens=uens;v4ens=uens;
% uens(:,1)=u(1)+0.1*randn(MM,1);
omega1=H1/sqrt(2);omega3=H2/sqrt(2);
    du=0.0125;dv1=0.0125;dv2=0.0125;dv3=0.0125;dv4=0.0125;
%     beta=2;
G1=zeros(5,5);%%% (u,v1,v2,v3,v4) to (u,pis1 pis-1,psi2,psi-2);
G1(1,1)=1;
G1(2,2)=1/2/sqrt(2)*(-1-1i);G1(2,3)=1/2/sqrt(2)*(1-1i);
G1(3,2)=1/2/sqrt(2)*(-1+1i);G1(3,3)=1/2/sqrt(2)*(1+1i);
G1(4,4)=1/2/sqrt(2)*(-1-1i);G1(4,5)=1/2/sqrt(2)*(1-1i);
G1(5,4)=1/2/sqrt(2)*(-1+1i);G1(5,5)=1/2/sqrt(2)*(1+1i);
% temp=inv(G1)\[u(:,1);psik(1,1);conj(psik(1,1));psik(2,1);conj(psik(2,1)) ];
    kk1dfull=[1,-1,2,-2;0,0,0,0];
    repk1d=repmat([1,-1,2,-2],L,1)*1i;
kk00=[];
temp=[u(i-1);psik(1,i-1);conj(psik(1,i-1));psik(2,i-1);conj(psik(2,i-1))];
for i=2:N
        if mod(i,floor(N/5)) == 0
%         disp(i*dt)
 fprintf('rest step is %d\n',N-i);
        end
 x_loc = [x(:,i-1),y(:,i-1)];%+sig_ex*sqrt(dt)*randn(L,2);  
Gy1d = [zeros(L,1),exp(1i * x_loc *  kk1dfull).*repk1d ];
Gu=[ones(L,1),zeros(L,4)];
Gux=Gu*G1;Guy=(Gu+Gy1d)*G1;
gx1=Gux*[uens(:,i-1),v1ens(:,i-1),v2ens(:,i-1),v3ens(:,i-1),v4ens(:,i-1)]';
gx2=Guy*[uens(:,i-1),v1ens(:,i-1),v2ens(:,i-1),v3ens(:,i-1),v4ens(:,i-1)]';
gx=[gx1;gx2];gx=gx';
if max(abs(imag(gx(:))))>10^(-10)
    disp('error complex vel')
end
gx=real(gx);
    temp = (gx - repmat(mean(gx),MM,1));
    PGu = temp'*((uens(:,i-1) - mean(uens(:,i-1))))/(MM-1)/(sig_ex^2);
    PGv1 = temp'*((v1ens(:,i-1) - mean(v1ens(:,i-1))))/(MM-1)/(sig_ex^2);
    PGv2 = temp'*((v2ens(:,i-1) - mean(v2ens(:,i-1))))/(MM-1)/(sig_ex^2);
    PGv3 = temp'*((v3ens(:,i-1) - mean(v3ens(:,i-1))))/(MM-1)/(sig_ex^2);
    PGv4 = temp'*((v4ens(:,i-1) - mean(v4ens(:,i-1))))/(MM-1)/(sig_ex^2);
 
% dz=[x(:,i);y(:,i)]-[x(:,i-1);y(:,i-1)];
    diff_x1 = x(:,i) - x(:,i-1); diff_x2 = x(:,i) - x(:,i-1) + L1; diff_x3 = x(:,i) - x(:,i-1) - L1;  
    diff_y1 = y(:,i) - y(:,i-1); diff_y2 = y(:,i) - y(:,i-1) + L1; diff_y3 = y(:,i) - y(:,i-1) - L1;  
    diff_xtemp = min(abs(diff_x1), abs(diff_x2)); diff_x_index = min(abs(diff_x3), diff_xtemp);
    diff_ytemp = min(abs(diff_y1), abs(diff_y2)); diff_y_index = min(abs(diff_y3), diff_ytemp);
    diff_x1_index = (diff_x_index == abs(diff_x1)); diff_x2_index = (diff_x_index == abs(diff_x2)); diff_x3_index = (diff_x_index == abs(diff_x3)); 
    diff_y1_index = (diff_y_index == abs(diff_y1)); diff_y2_index = (diff_y_index == abs(diff_y2)); diff_y3_index = (diff_y_index == abs(diff_y3)); 
    diff_x = diff_x1 .* diff_x1_index + diff_x2 .* diff_x2_index + diff_x3 .* diff_x3_index;
    diff_y = diff_y1 .* diff_y1_index + diff_y2 .* diff_y2_index + diff_y3 .* diff_y3_index;
    dz = [diff_x; diff_y];
dz=(dz*ones(1,MM))';

DA_term_u=-PGu'*(gx*dt-dz+sig_ex*sqrt(dt)*randn(MM,2*L))';
DA_term_v1=-PGv1'*(gx*dt-dz+sig_ex*sqrt(dt)*randn(MM,2*L))';
DA_term_v2=-PGv2'*(gx*dt-dz+sig_ex*sqrt(dt)*randn(MM,2*L))';
DA_term_v3=-PGv3'*(gx*dt-dz+sig_ex*sqrt(dt)*randn(MM,2*L))';
DA_term_v4=-PGv4'*(gx*dt-dz+sig_ex*sqrt(dt)*randn(MM,2*L))';
% gbar=repmat(mean(gx),MM,1);
% DA_term_u=-PGu'/2*(gx*dt+gbar*dt-2*dz)';
% DA_term_v1=-PGv1'/2*(gx*dt+gbar*dt-2*dz)';
% DA_term_v2=-PGv2'/2*(gx*dt+gbar*dt-2*dz)';
% DA_term_v3=-PGv3'/2*(gx*dt+gbar*dt-2*dz)';
% DA_term_v4=-PGv4'/2*(gx*dt+gbar*dt-2*dz)';

uens(:,i)=uens(:,i-1)+(omega1*v1ens(:,i-1)+2*omega3*v3ens(:,i-1)-du*uens(:,i-1))*dt+...
    sigmau*randn(MM,1)*sqrt(dt)+DA_term_u';
v1ens(:,i)=v1ens(:,i-1)+(-beta*v2ens(:,i-1)+v2ens(:,i-1).*uens(:,i-1)-...
    2*omega1*uens(:,i-1)-dv1*v1ens(:,i-1))*dt+sigmau*randn(MM,1)*sqrt(dt)+DA_term_v1';
v2ens(:,i)=v2ens(:,i-1)+(beta*v1ens(:,i-1)+v1ens(:,i-1).*uens(:,i-1)-...
    dv2*v2ens(:,i-1))*dt+sigmau*randn(MM,1)*sqrt(dt)+DA_term_v2';  
v3ens(:,i)=v3ens(:,i-1)+(-beta/2*v4ens(:,i-1)+2*v4ens(:,i-1).*uens(:,i-1)-...
     omega3*uens(:,i-1)-dv3*v3ens(:,i-1))*dt+sigmau*randn(MM,1)*sqrt(dt)+DA_term_v3'; 
v4ens(:,i)=v4ens(:,i-1)+(beta/2*v3ens(:,i-1)-2*v3ens(:,i-1).*uens(:,i-1)-...
    dv4*v4ens(:,i-1))*dt+sigmau*randn(MM,1)*sqrt(dt)+DA_term_v4';   
end

psi1ens=1/2/sqrt(2)*( (v2ens-v1ens)-1i*(v2ens+v1ens) );
psi2ens=1/2/sqrt(2)*( (v4ens-v3ens)-1i*(v4ens+v3ens) );
psi1=mean(psi1ens);psi2=mean(psi2ens);
uass=mean(uens);
estuhat2=zeros(1+2*Kmax,N);
estuhat2(1,:)=uass;
estuhat2(2,:)=psi1;estuhat2(3,:)=conj(psi1);
estuhat2(4,:)=psi2;estuhat2(5,:)=conj(psi2);

uexact=u(:,1/dt:1/dt:end);
psikexact=psik(:,1/dt:1/dt:end);
estuhat=estuhat2(:,1/dt:1/dt:end);

uhate=zeros(1+2*Kmax,N);
uhate(1,:)=u;
uhate(2:2:end,:)=psik;
uhate(3:2:end,:)=conj(psik);


% rmscc(estuhat2(1,:),uhate(1,:),1);
% rmscc(estuhat2(2,:),uhate(2,:),1);
% rmscc(estuhat2(4,:),uhate(4,:),1);
rec=1/dt:1/dt:N;
rmscc(estuhat2(1,rec),uhate(1,rec),1);
rmscc(estuhat2(2,rec),uhate(2,rec),1);
rmscc(estuhat2(4,rec),uhate(4,rec),1);

a=estuhat2(1,1/dt:1/dt:end);save a a
b=uhate(1,1/dt:1/dt:end);save b b
uhate=uhate(:,1/dt:1/dt:end);

compute_errorlayer1d
  errest=[mean(rmstime(10:end)),mean(rmsall),mean(ccall)];
clear x y u  psik 
% save uhate uhate;save estuhat2 estuhat2;
fprintf('test case is %d\n',testcase);
% uass=uass(1:1000:end);
% uex=u(1:1000:end);