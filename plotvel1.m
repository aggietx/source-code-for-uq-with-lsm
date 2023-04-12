Dim_Grid = 30;
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
    Ga = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(1,:))); % Fourier bases for u
    Gb = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(2,:))); % Fourier bases for v
datacase=2;ii=6;
if datacase==1
data=gamma_mean_trace;
% data=gamma_mean_trace_local;
% data=all_filmean(:,:,ii);
elseif datacase==2
data=gamma_mean_trace_smoother;
%  data=all_smmean(:,:,ii);
%  data=u_hat;
else
data=gamma_mean_trace_smoother_sample;
%  u_hatapp=all_bsmean(:,:,ii);
end
dx=L1/Dim_Grid;
%% exact
time=15;
    u = (Ga*  data(:,time/dt)); 
    v = (Gb*  data(:,time/dt)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity, record exact'  )
      max(abs(imag(u)))
    end
    
    u=real(u);v=real(v);
 u=reshape(u, Dim_Grid,Dim_Grid);
 v=reshape(v, Dim_Grid,Dim_Grid);
figure();
    imagesc([0,L1],[0,L1],sqrt(u.^2+v.^2));colorbar
%     temp1=caxis;
  mymap = 1-[0 0 0
1/10 0 0
0 1/10 0
0 0 1/10
1/10 1/10 1/10];
colormap(mymap)
set(gca,'ydir','normal')
% colormap(cmocean('balance'))

  hold on
  quiver(dx/2:dx:L1-dx/2,(dx/2:dx:L1-dx/2)',u,v,1.3,'color',[1 1 1]*0.4);
        set(gca,'fontsize',18)
%     xlabel('Km','FontSize',16)  ;ylabel('Km','FontSize',16)      
     
title(['Day ', num2str(time)],'FontSize',16) 