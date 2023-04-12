%%% compute exact layer vel and estimated vel
% get psikexact.mat  uexact.mat estuhat.mat 
%load psikexact; load uexact;load estuhat.mat 
%%%% 2d data
close all;
% load data72
% load data36;
% load data36_p05
% load data72_p05
load data
% load ccall;
K_max=(sqrt(size(estuhat,1))-1)/2;Kuse=K_max;if K_max~=Kuse;reduced=1;else reduced=0;end;
if Kuse>K_max;Kuse=K_max;reduced=0;end
dimuhat=(2*K_max+1)^2;Kmax=K_max;
[ky,kx]=meshgrid([0:K_max,-K_max:-1],[0:K_max,-K_max:-1]);
kk=[kx(:),ky(:)]';
rk=zeros(size(kk));
for i=2:dimuhat
        ki=kk(:,i);
%         rk(:,i)=1i*[-ki(2);ki(1)]/norm(ki);
        rk(:,i)=1i*[-ki(2);ki(1)]/sqrt(ki(1)^2+ki(2)^2+1);
end
% ix=min(K_max,k_0);iy=min(K_max,k_0);ind=getind(ix,iy,kk);%% sample plot
kk0=1:10;kk0=kk0';
redind=getreducedind(kk,K_max,Kuse);

S=3;%%% quiver scaling
nx=20;ny=nx;Dim_Grid = nx;L1=2*pi;

% ind=find(ccall<0.5);
% ind=find(ccall==max(ccall));
% ind=find(ccall>0.9);
% time=ind(1);
time=392;
x1d=linspace(0,L1,nx);x1d=x1d';
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
    Ga = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(1,:))); % Fourier bases for u
    Gb = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(2,:))); % Fourier bases for v
    uxest = (Ga* estuhat(:,time)); 
    uyest = (Gb* estuhat(:,time)); 
    if max(abs(imag(uxest)))>10^(-6) || max(abs(imag(uyest)))>10^(-6)
      disp('complex velocity, record exact'  )
      max(abs(imag(uxest)))
    end
    
    uxest=real(uxest);uyest=real(uyest);
 uxest=reshape(uxest, Dim_Grid,Dim_Grid);
 uyest=reshape(uyest, Dim_Grid,Dim_Grid);
 figure()
 subplot(2,1,2)
     imagesc([0,L1],[0,L1],sqrt(uxest.^2+uyest.^2));colorbar
%     temp1=caxis;
lightcolor
hold on

    quiver(xx, yy, uxest, uyest,S, 'linewidth',1.5)
    xlim([0, L1 ])
    ylim([0, L1 ])  
    xlabel('x','FontSize',18);  ylabel('y','FontSize',18)
    setgca(18)
        title(['Estimated velocity at t = ', num2str(time)],'FontSize',24)


for i=time; 
a=psikexact(:,i);b=real(a)-1i*imag(a);
a1=1i.*a.*kk0;b1=-1i*b.*kk0;
uy1d = exp(1i * x1d * kk0')*a1+exp(-1i * x1d * kk0')*b1;
uy=repmat(uy1d',ny,1)+uexact(i);
ux=ones(ny,nx)*uexact(i);
 subplot(2,1,1)

     imagesc([0,L1],[0,L1],sqrt(ux.^2+uy.^2));colorbar
%     temp1=caxis;
lightcolor
hold on
    quiver(xx, yy, ux, uy, S,'linewidth',1.5)
%     hold on
%     scatter(x(:,i),y(:,i),'red','linewidth',2);
    xlim([0, L1 ])
    ylim([0, L1 ])  
    xlabel('x','FontSize',18);  ylabel('y','FontSize',18)
    setgca(18)
    box on    
    title(['Exact velocity field t = ', num2str(i)],'FontSize',24)

end
fprintf('time is  %2.1f\n',time);
rmscc(sqrt(uxest.^2+uyest.^2),sqrt(ux.^2+uy.^2),1);
fprintf('Kmax is  %d\n',K_max);
