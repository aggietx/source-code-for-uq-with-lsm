clear;
% close all;
LL=[6,12:12:96];
ent1=[34.60 23.30
     24.63 38.33
     15.77 58.83
     11.59 72.62
     9.24 82.71
     7.77 90.43
     6.79 96.78
     6.05 102.06
     5.53 106.52];err1=[0.64  0.77 
                       0.52  0.85
                       0.40  0.92
                       0.34  0.94
                       0.30  0.95
                        0.27  0.96
                        0.25  0.97
                        0.24  0.97
                        0.23  0.97];
ent2=[37.91 29.33
     26.49 42.13
    16.56 59.98
    12.05 72.04
    9.59 81.54
    8.06 89.05
    7.04 95.36
    6.25 100.53
    5.72 104.88            ];
                    err2=[ 0.67  0.75 
                        0.52  0.85
                        0.41  0.91
                        0.34  0.94
                        0.30  0.95
                        0.28  0.96
                        0.26  0.97
                        0.24  0.97
                        0.23,0.97];
   err3=[0.40  0.92
       0.33  0.95
       0.30  0.96
        0.30  0.96
        0.30  0.96
        0.30  0.96
         0.30  0.96
         0.29  0.96
         0.29  0.96];% close all;
  
     figure();
     subplot(2,2,1);
     
     plot(LL,ent1(:,1)','-*k','linewidth',2);hold on
    plot(LL,ent2(:,1)','-*b','linewidth',2);hold on
    title('(a) Signal','fontsize',16);
    xlim([6,96]);setgca(16);
    
    set(gca,'XTick',LL,'fontsize',16)
    
      subplot(2,2,2);
           plot(LL,ent1(:,2)','-*k','linewidth',2);hold on
    plot(LL,ent2(:,2)','-*b','linewidth',2);hold on
    title('(b) Dispersion','fontsize',16);
    xlim([6,96]);setgca(16);
       set(gca,'XTick',LL,'fontsize',16)
    
         subplot(2,2,3);
     
     plot(LL,err1(:,1)','-*k','linewidth',2);hold on
    plot(LL,err2(:,1)','-*b','linewidth',2);hold on
     plot(LL,err3(:,1)','-*r','linewidth',2);hold on
    title('(c) RMS','fontsize',16);
    xlim([6,96])
    setgca(16);
    xlabel('L','fontsize',16);
       set(gca,'XTick',LL,'fontsize',16)
       
      subplot(2,2,4);
           plot(LL,err1(:,2)','-*k','linewidth',2);hold on
    plot(LL,err2(:,2)','-*b','linewidth',2);hold on
      plot(LL,err3(:,2)','-*r','linewidth',2);hold on
    title('(d) CC','fontsize',16);
    xlim([6,96]);setgca(16);
       set(gca,'XTick',LL,'fontsize',16)
    lgnd=legend('Exact parameters','Estimated parameters','Free run with fixed noise');
set(lgnd,'FontSize',16);
xlabel('L','fontsize',16);