Dim_Grid = 48;
[xx,yy] = meshgrid(linspace(0,L1,Dim_Grid), linspace(0,L1,Dim_Grid));
xy=[xx(:),yy(:)];
    Ga = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(1,:))); % Fourier bases for u
    Gb = (exp(1i * xy * kk*2*pi/L1) .* (ones(Dim_Grid^2,1) * rk(2,:))); % Fourier bases for v

dx=L1/Dim_Grid;
%% exact
time=15;
    u = (Ga*  u_hat(:,time/dt)); 
    v = (Gb*  u_hat(:,time/dt)); 
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
     
title(['Exact flow at day ', num2str(time)],'FontSize',16) 

%% filter
    u = (Ga*  gamma_mean_trace(:,time/dt)); 
    v = (Gb*  gamma_mean_trace(:,time/dt)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity, record filter'  )
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
     
title(['Filtered flow at day ', num2str(time)],'FontSize',16) 
if exist('gamma_mean_trace_smoother','var')
%% smoother
    u = (Ga*  gamma_mean_trace_smoother(:,time/dt)); 
    v = (Gb*  gamma_mean_trace_smoother(:,time/dt)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity, record smoother'  )
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
     
title(['Flow with smoother at day ', num2str(time)],'FontSize',16) 


end


%% backward sampling
if exist('gamma_mean_trace_smoother_sample','var')
    u = (Ga*  gamma_mean_trace_smoother_sample(:,time/dt)); 
    v = (Gb*  gamma_mean_trace_smoother_sample(:,time/dt)); 
    if max(abs(imag(u)))>10^(-6) || max(abs(imag(v)))>10^(-6)
      disp('complex velocity, record backward sampling'  )
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
     
title(['Flow with backward sampling at day ', num2str(time)],'FontSize',16) 

end

