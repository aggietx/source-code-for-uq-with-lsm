clear ;
% close all;
%%%%%% should compare diagcov, constant cov and full cov (memory huge)
%%%% 36,72,108 
p=.5;
if p==1
rmsL=[6.32,2.24,1.22,0.88,0.76,0.75,0.70,0.69,0.68,0.68; 
    3.7,1.14,0.72,0.62,0.63,0.61,0.60,0.61,0.60,0.61;
    2.27,0.84,0.65,0.70,0.58,0.56,0.56,0.56,0.57,0.56];
ccL=[0.06,0.20,0.44,0.63,0.72,0.73,0.77,0.78,0.78,0.78;
   0.15,0.49,0.75,0.82,0.82, 0.83,0.83,0.83,0.83,0.83;
   0.28,0.67,0.79,0.77,0.85,0.86,0.86,0.86,0.86,0.86];
else
 rmsL=[2.40,0.96,0.63,0.53,0.50,0.49,0.50,0.49,0.49,0.50;
     1.27,0.55,0.40,0.38, 0.37, 0.37,0.37,0.37, 0.37,0.37 ;
     1.05,0.49,0.41,0.39 ,0.38,0.38,0.38,0.38,0.38,0.38];
 ccL=[0.21,0.59,0.79,0.86,0.87,0.88,0.88,0.88,0.88,0.88;
     0.47,0.84,0.91,0.93,0.93,0.93,0.93,0.93,0.93,0.93;
      0.58, 0.87,0.91,0.92,0.92,0.93,0.92,0.93,0.93,0.93];
end
figure();
subplot(2,1,1)
plot(1:10,rmsL(1,:),'-*r','linewidth',2);hold on;
plot(1:10,rmsL(2,:),'-*g','linewidth',2);hold on;
plot(1:10,rmsL(3,:),'-*b','linewidth',2);hold on;
title('Average RMS of the velocity','fontsize',18)
setgca(18)
xlim([1,10]);
subplot(2,1,2)
plot(1:10,ccL(1,:),'-*r','linewidth',2);hold on;
plot(1:10,ccL(2,:),'-*g','linewidth',2);hold on;
plot(1:10,ccL(3,:),'-*b','linewidth',2);hold on;
xlim([1,10])
title('Average CC of the velocity','fontsize',18)
setgca(18)
xlabel('Iteration','fontsize',24)
lgnd=legend('L=36','L=72','L=108');
set(lgnd,'FontSize',18);
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load data36_p05;
figure();

subplot(2,1,1);

 plot(abs(uexact),'*-','linewidth',1)
 title('Absolute uexact','FontSize',24)
setgca(18);
subplot(2,1,2);
plot(ccall,'*-','linewidth',1)
setgca(18);
title('cc','FontSize',24)
xlabel('t','fontsize',24)
xlim([1,length(uexact)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare using different Kmax
%%% L=36, Kmax=10 or 12 or 6;
if p==1
rmsK=[2.19,1.21,1.14 ,1.13,1.13,1.13,1.13,1.13,1.13,1.13;
    6.32,2.24,1.22,0.88,0.76,0.75,0.70,0.69,0.68,0.68;  %%%6.55,2.10,1.11,0.87,0.8,0.77,0.75,0.73,0.71,0.70;%%%T=200 
   6.79,2.43,1.31,0.92,0.79,0.73,0.74,0.73,0.73,0.73];
ccK=[0.10,0.25,0.32,0.34, 0.35,0.35,0.35,0.36,0.36,0.36;
    0.06,0.20,0.44,0.63,0.72,0.73,0.77,0.78,0.78,0.78;%%%0.04,0.19,0.44,0.62,0.69,0.71,0.73,0.74,0.75,0.76;
    0.04,0.16,0.38,0.59,0.70,0.74,0.74,0.74,0.74,0.74];
else
rmsK=[1.47,1.26,1.26,1.26,1.26,1.26,1.27,1.27,1.27,1.27;
    2.40,0.96,0.63,0.53,0.50,0.49,0.50,0.49,0.49,0.50;
    2.82,1.11,0.68,0.55,0.51,0.50,0.51,0.50,0.51,0.50];
ccK=[0.10,0.19,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2;
    0.21,0.59,0.79,0.86,0.87,0.88,0.88,0.88,0.88,0.88;
    0.17,0.52,0.76,0.84,0.87,0.87,0.87,0.87,0.87,0.87];
end
figure();
subplot(2,1,1)
plot(1:10,rmsK(1,:),'-*r','linewidth',2);hold on;
plot(1:10,rmsK(2,:),'-*g','linewidth',2);hold on;
plot(1:10,rmsK(3,:),'-*b','linewidth',2);hold on;
title('Average RMS of the velocity','fontsize',18)
setgca(18)
xlim([1,10]);
subplot(2,1,2)
plot(1:10,ccK(1,:),'-*r','linewidth',2);hold on;
plot(1:10,ccK(2,:),'-*g','linewidth',2);hold on;
plot(1:10,ccK(3,:),'-*b','linewidth',2);hold on;
xlim([1,10])
title('Average CC of the velocity','fontsize',18)
setgca(18)
xlabel('Iteration','fontsize',24)
lgnd=legend('Kmax=6','Kmax=10','Kmax=12');
set(lgnd,'FontSize',18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%