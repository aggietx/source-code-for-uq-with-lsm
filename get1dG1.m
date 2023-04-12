function [G1,G2]=get1dG1(Kmax,x_loc )
%%%% for layered model compute vel
kk1dfull=zeros(2,2*Kmax);
for i=1:Kmax
   kk1dfull(1,2*i-1)=i;
   kk1dfull(1,2*i)=-i;
end
% [xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
% x_loc=[xx(:),yy(:)];
dim=size(x_loc,1);
repk1d2=repmat(kk1dfull(1,:),dim,1)*1i;
Gy1d = [zeros(dim,1),exp(1i * x_loc *  kk1dfull).*repk1d2 ];
Gu=[ones(dim,1),zeros(dim,2*Kmax)];
G1=Gu;G2=(Gu+Gy1d);
