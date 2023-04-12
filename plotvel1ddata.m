%%%% plot 1d vel comparison
% save estuhat estuhat;save uhate uhate
% get estuhat.mat uhate.mat
%%% plot(
close all;
% load data1d_p05;
load data1d_p1
load xe;load ye;

type=3;%%% u: 1:large: x flow 2: medium 3,small block
num=1;%%% which one 
% num=3;
S1=30;u0=real(uexact(1,S1:end));
u=real(uexact(1,:));
if type==1
ind=find(u0>.6);
% ind=find(u0>1);
elseif type==2
         ind=find(abs(u0-0.3)<0.05);
%     ind=find(abs(u0-0.5)<0.05);
else
ind=find(abs(u0)<0.05);
end
ind=ind+S1-1;
time=ind(num);
fprintf('type %d, total time %d, mean rms and mean cc %2.3f %2.3f\n',type,length(ind),mean(rmsall(ind)),mean(ccall(ind)));

% time=392;
Kmax=(size(estuhat,1)-1)/2;S=3;
% close all;
nx=20;ny=nx;Dim_Grid = nx;L1=2*pi;
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
x_loc=[xx(:),yy(:)];
kk1dfull=zeros(2,2*Kmax);
for i=1:Kmax
   kk1dfull(1,2*i-1)=i;
   kk1dfull(1,2*i)=-i;
end
kk1dfullf=[[0;0],kk1dfull];
repk1d=repmat(kk1dfull(1,:),Dim_Grid^2,1)*1i;

Gy1d = [zeros(Dim_Grid^2,1),exp(1i * x_loc *  kk1dfull).*repk1d ];
Gu=[ones(Dim_Grid^2,1),zeros(Dim_Grid^2,2*Kmax)];
G1=Gu;G2=(Gu+Gy1d);   
velx=G1*uhate(:,time);vely=G2*uhate(:,time);
velx=reshape(velx,nx,ny);vely=reshape(vely,nx,ny);
velxa=G1*estuhat(:,time);velya=G2*estuhat(:,time);
velxa=reshape(velxa,nx,ny);velya=reshape(velya,nx,ny);
figure();
subplot(2,1,1);imagesc([0,L1],[0,L1],sqrt(velx.^2+vely.^2));colorbar; lightcolor; hold on
% quiver(velx,vely)
    quiver(xx, yy, velx,vely, S,'linewidth',1.5)
            hold on
    scatter(xe(:,time),ye(:,time),'red','linewidth',2);
    xlabel('x','FontSize',18);  ylabel('y','FontSize',18)

    xlim([0, L1 ])
    ylim([0, L1 ])  
     setgca(18)
             title(['Exact velocity at t = ', num2str(time)],'FontSize',24)

subplot(2,1,2);
imagesc([0,L1],[0,L1],sqrt(real(velxa).^2+real(velya).^2));colorbar; lightcolor;hold on
% quiver(velxa,velya);
    quiver(xx, yy, velxa,velya, S,'linewidth',1.5)
            hold on
    scatter(xe(:,time),ye(:,time),'red','linewidth',2);
    xlabel('x','FontSize',18);  ylabel('y','FontSize',18)

rmscc(sqrt(real(velxa).^2+real(velya).^2),sqrt(velx.^2+vely.^2),1);
    xlim([0, L1 ])
    ylim([0, L1 ])  
    setgca(18)
            title(['Estimated velocity at t = ', num2str(time)],'FontSize',24)

% figure
% subplot(2,1,1);
% plot(uhate(1,:));
% subplot(2,1,2);
% plot(real(estuhat(1,:)));

% plot2line(uhate(1,:),real(estuhat(1,:)));