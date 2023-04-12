    
function [uv,vu]=dyad_infolocal(uinit,vinit, Ti, MM,c,F_u,sigma_u,d_v,sigma_v,dt,d_u)
uv=zeros(Ti/dt,1);
vu=zeros(Ti/dt,1);
uold=uinit*ones(MM,1);
vold=vinit*ones(MM,1);
    for i=2:Ti/dt
        v=vold;u=uold;
    unew=uold+dt*((-d_u+c*v).*u+F_u)+sigma_u*sqrt(dt)*randn(MM,1);
    vnew=vold+(-d_v*v-c*u.^2)*dt+sigma_v*sqrt(dt)*randn(MM,1);   
    mu1=mean(unew);mu2=mean(vnew);
diffx1=unew-mu1;diffx2=vnew-mu2;
r11=mean(diffx1.*diffx1);r12=mean(diffx1.*diffx2);
r21=r12;r22=mean(diffx2.*diffx2);  
        vu(i)=2*c*mu1*r21/r22;
        uv(i)=-2*c*mu2-c*mu1*r12/r11;
        vold=vnew;
        uold=unew;
    end