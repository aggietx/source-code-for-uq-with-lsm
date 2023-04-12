%%%% plot vel for topography

nx=20;ny=nx;Dim_Grid = nx;L1=2*pi;
x=linspace(0,L1,nx);x=x';
time=756120;
for i=time; 
a=psik(:,i);b=real(a)-1i*imag(a);
a1=1i.*a.*kk0;b1=-1i*b.*kk0;
uy1d = exp(1i * x * kk0')*a1+exp(-1i * x * kk0')*b1;
end
uy=repmat(uy1d',ny,1)+u(i);
ux=ones(ny,nx)*u(i);

save ux ux
save uy uy

% load ux;load uy;
close  all;
figure();quiver(ux,uy)

return
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];kk=[0;0];
    Ga = (exp(1i * xy * kk*2*pi/L1)); % Fourier bases for u
    Gb = (exp(1i * xy * kk*2*pi/L1)); % Fourier bases for v
ux=Ga*u(i);uy=Gb*u(i);
ux=reshape(ux,Dim_Grid,Dim_Grid);
uy0=reshape(uy,Dim_Grid,Dim_Grid);
uy=uy0+uy2d;

