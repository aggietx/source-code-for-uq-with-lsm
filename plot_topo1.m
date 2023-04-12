close all;clear
load acfpsi1;load  acfpsi75;load  acfpsi5;
load psi1;load psi75;load psi5;
load streamfun1;load streamfun5;load streamfun75;
psi1=real(psi1);psi5=real(psi5);psi75=real(psi75);
time1=200;time2=300;
plotcase=3;
if plotcase==1
psi=psi1;acfpsi=acfpsi1;streamfun=streamfun1;
elseif plotcase==2
psi=psi75;acfpsi=acfpsi75;streamfun=streamfun75;
elseif plotcase==3
psi=psi5;acfpsi=acfpsi5;streamfun=streamfun5;
end
ts=150;te=350;dt=1/1000;gap=20;dt0=gap*dt;
figure();
subplot(7,4,1);
plot(ts:dt0:te,psi(end,ts/dt0:te/dt0),'b','linewidth',2)
xlim([ts,te]);setgca(14);
ylabel('u','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
title('(a) Trajectory','fontsize',14)
Kmax=10;
for i=1:Kmax
    j0=floor(i/2);
    if i==1;
        j=5;
    elseif i==2
        j=6;
    elseif i==3
        j=9;
    elseif i==4
        j=10;    
    elseif i==5
        j=13;
    elseif i==6
        j=14;
    elseif i==7
        j=17; 
    elseif i==8
        j=18;
    elseif i==9
        j=21;
    elseif i==10
        j=22;         
    end
  subplot(7,4,j);
plot(ts:dt0:te,psi(i,ts/dt0:te/dt0),'b','linewidth',2)
xlim([ts,te]);setgca(14);

    if i==1;
ylabel('\psi_1','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==2
ylabel('\psi_2','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==3
ylabel('\psi_3','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==4
ylabel('\psi_4','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
  
    elseif i==5
ylabel('\psi_5','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==6
ylabel('\psi_6','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==7
ylabel('\psi_7','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==8
ylabel('\psi_8','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==9
ylabel('\psi_9','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==10
ylabel('\psi_{10}','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
       
    end
    
    subplot(7,4,j+2);
plot(0:dt:60,acfpsi(i,:),'b','linewidth',2)
xlim([0,60]);setgca(14);
    if i==1;
ylabel('\psi_1','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==2
ylabel('\psi_2','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==3
ylabel('\psi_3','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==4
ylabel('\psi_4','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
  
    elseif i==5
ylabel('\psi_5','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==6
ylabel('\psi_6','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==7
ylabel('\psi_7','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==8
ylabel('\psi_8','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==9
ylabel('\psi_9','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)

    elseif i==10
ylabel('\psi_{10}','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
       
    end
end
subplot(7,4,3);
plot(0:dt:60,acfpsi(end,:),'b','linewidth',2)
xlim([0,60]);setgca(14);
ylabel('u','fontsize',14);set(get(gca,'YLabel'),'Rotation',0)
title('(b) ACF','fontsize',14)
return
nx=100;
x=0:2*pi/nx:2*pi;x=x';

subplot(7,4,25);
plot(x,streamfun(:,time1),'r','linewidth',2)
xlim([0,2*pi]);setgca(14);
title(['(c) stream function at t = ', num2str(time1)],'FontSize',14) 

subplot(7,4,27);
plot(x,streamfun(:,time2),'r','linewidth',2)
xlim([0,2*pi]);setgca(14);
title(['(e) stream function at t = ', num2str(time2)],'FontSize',14) 
