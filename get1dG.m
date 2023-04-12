function [G1,G2]=get1dG(Kmax,Dim_Grid,L1)
%%%% for layered model compute vel
kk1dfull=zeros(2,2*Kmax);
for i=1:Kmax
   kk1dfull(1,2*i-1)=i;
   kk1dfull(1,2*i)=-i;
end
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
x_loc=[xx(:),yy(:)];

repk1d2=repmat(kk1dfull(1,:),Dim_Grid^2,1)*1i;
Gy1d = [zeros(Dim_Grid^2,1),exp(1i * x_loc *  kk1dfull).*repk1d2 ];
Gu=[ones(Dim_Grid^2,1),zeros(Dim_Grid^2,2*Kmax)];
G1=Gu;G2=(Gu+Gy1d);
