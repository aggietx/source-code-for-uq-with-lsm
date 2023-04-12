
close all;clear
p=1;L=36;maxu=0.4;
type=1;num=1;S=2;%%% quiver scaling

if L==36 
    if p==0.5
load xye36_p05;load data_p05_L36_para2d;
    else
load xye36_p1;load data_p1_L36_para2d;        
    end
else
    if p==0.5
     load xye12_p05;load data_p05_L12_para2d;
    else
        load xye12_p1;load data_p1_L12_para2d; 
    end
end
xe=xye(1:L,:);
ye=xye(1+L:end,:);

S1=30;u0=real(uexact(1,S1:end));
u=real(uexact(1,:));
if type==1
ind=find(abs(u0)>maxu);
% ind=find(u0>1);
elseif type==2
         ind=find(abs(abs(u0)-maxu/2)<0.05);
%     ind=find(abs(u0-0.5)<0.05);
else
ind=find(abs(u0)<0.05);
end
ind=ind+S1-1;
time=ind(num);time1=time;

[ux,uy]=compute1dvel(Kmax,L1,Dim_Grid,uhate(:,time));
%% exact vel
 figure();
subplot(3,1,1);imagesc([0,L1],[0,L1],sqrt(ux.^2+uy.^2));colorbar; lightcolor; hold on
% quiver(velx,vely)
    quiver(xx, yy, ux,uy, S,'linewidth',1.5)
            hold on
    scatter(xe(:,time),ye(:,time),'red','linewidth',2);
%     xlabel('x','FontSize',14);  ylabel('y','FontSize',16)

    xlim([0, L1 ])
    ylim([0, L1 ])  
     setgca(14)
title(['Exact velocity  at t = ', num2str(time)],'FontSize',16)

     %% 2d estimation
     [ux2d,uy2d]=computeGBvel(rk,kk,L1,Dim_Grid,estuhat(:,time));

     
subplot(3,1,3);imagesc([0,L1],[0,L1],sqrt(ux2d.^2+uy2d.^2));colorbar; lightcolor; hold on
% quiver(velx,vely)
    quiver(xx, yy, ux2d,uy2d, S,'linewidth',1.5)
            hold on
    scatter(xe(:,time),ye(:,time),'red','linewidth',2);
    xlabel('x','FontSize',14);  ylabel('y','FontSize',16)

    xlim([0, L1 ])
    ylim([0, L1 ])  
     setgca(14)
             title(['Estimated velocity (2D) at t = ', num2str(time)],'FontSize',16)             

% fprintf('2d type %d, total time %d, mean rms and mean cc %2.3f %2.3f\n',type,length(ind),mean(rmsall(ind)),mean(ccall(ind)));
fprintf('2d type %d, total time %d, mean rms %2.3f\n',type,length(ind),mean(rmstime(ind)));
     %% 1d estimation
if L==36 
    if p==0.5
load data_p05_L36_para1d;
    else
load data_p1_L36_para1d;        
    end
else
    if p==0.5
load data_p05_L12_para1d;
    else
load data_p1_L12_para1d; 
    end
end

time=time1;
 [ux1d,uy1d]=compute1dvel(Kmax,L1,Dim_Grid,estuhat(:,time));
%  fprintf('1d type %d, total time %d, mean rms and mean cc %2.3f %2.3f\n',type,length(ind),mean(rmsall(ind)),mean(ccall(ind)));
fprintf('1d type %d, total time %d, mean rms %2.3f\n',type,length(ind),mean(rmstime(ind)));


subplot(3,1,2);imagesc([0,L1],[0,L1],sqrt(ux1d.^2+uy1d.^2));colorbar; lightcolor; hold on
% quiver(velx,vely)
    quiver(xx, yy, ux1d,uy1d, S,'linewidth',1.5)
            hold on
    scatter(xe(:,time),ye(:,time),'red','linewidth',2);
%     xlabel('x','FontSize',14);  ylabel('y','FontSize',16)

    xlim([0, L1 ])
    ylim([0, L1 ])  
     setgca(14)
             title(['Estimated velocity (1D) at t = ', num2str(time)],'FontSize',16)
