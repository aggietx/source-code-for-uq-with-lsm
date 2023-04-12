function [velx,vely,xx,yy]=compute1dvel(Kmax,L1,Dim_Grid,uhat)
% [velx,vely]=compute1dvel(Kmax,L1,Dim_Grid,estuhat(:,100))
% nx=20;ny=nx;Dim_Grid = nx;L1=2*pi;
nx=Dim_Grid;ny=nx;
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
velx=G1*uhat;vely=G2*uhat;
velx=reshape(velx,nx,ny);vely=reshape(vely,nx,ny);
% whos velx
    if max(abs(imag(velx(:))))>10^(-6) || max(abs(imag(vely(:))))>10^(-6)
      disp('complex velocity, record exact'  )
      max(abs(imag(velx)))
    end
   velx=real(velx);
vely=real(vely);
