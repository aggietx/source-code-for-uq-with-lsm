function [ux,uy,xx,yy]=computeGBvel(rk,kk,L1,Dim_Grid,uhat)
% % [ux,uy]=computeGBvel(rk,kk,L1,20,u_hat(:,100))

[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
    Ga = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(1,:))); % Fourier bases for u
    Gb = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(2,:))); % Fourier bases for v
    ux = (Ga* uhat); 
    uy = (Gb* uhat); 
    if max(abs(imag(ux)))>10^(-6) || max(abs(imag(uy)))>10^(-6)
      disp('complex velocity'  )
      max(abs(imag(ux)))
    end
    
    ux=real(ux);uy=real(uy);
 ux=reshape(ux, Dim_Grid,Dim_Grid);
 uy=reshape(uy, Dim_Grid,Dim_Grid);
