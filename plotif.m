%%%% plot information flow and mean
% close all;
figure
subplot(3,2,1);
% plot(0:dt:T,Tf,'r',...
%     0:dt:T,Tn,'--b','linewidth',2)
% lgnd=legend('T^f'  ,  'T^n');
plot(0:dt:T,Tf,'r','linewidth',2)
lgnd=legend('T^f' );
set(lgnd,'FontSize',14,'Location', 'Best');
set(gca,'FontSize',16)
xlabel('time','FontSize',16) 
ylabel('info flow','FontSize',16) 

if testcase==1
% title(['stable node,  ', '\alpha is ', num2str(alpha), ',  \beta is ', num2str(beta)],'FontSize',16) 
title(['\alpha= ', num2str(alpha), ',  \beta= ', num2str(beta), ',  \sigma_{22}= ', num2str(sigma22)],'FontSize',16) 

%%% (neg,neg):stable; (neg,pos):saddle
elseif testcase==2
% title('unstable node','FontSize',16) 
title(['unstable node,  ', '\alpha is ', num2str(alpha), ',  \beta is ', num2str(beta)],'FontSize',16) 
%%%  (pos,pos):unstable;
elseif testcase==3
% title('saddle node','FontSize',16) 
title(['saddle node,  ', '\alpha is ', num2str(alpha), ',  \beta is ', num2str(beta)],'FontSize',16) 

%%% (neg,pos):saddle;
elseif testcase==4
% title('center','FontSize',16) 
title(['center,  ', '\alpha is ', num2str(alpha), ',  \beta is ', num2str(beta)],'FontSize',16) 

%%% (neg,pos):center
elseif testcase==5
% title('stable spiral','FontSize',16) 
title(['stable spiral,  ', '\alpha is ', num2str(alpha), ',  \beta is ', num2str(beta)],'FontSize',16) 

%%% neg+a*sqrt(-1):stable spiral
elseif testcase==6  
% title('unstable spiral','FontSize',16) %%% neg+a*sqrt(-1):unstable spiral   
title(['unstable spiral,  ', '\alpha is ', num2str(alpha), ',  \beta is ', num2str(beta)],'FontSize',16) 

end

subplot(3,2,2);
% plot(0:dt:T,Tf,'r',...
%     0:dt:T,Tn,'--b','linewidth',2)
% lgnd=legend('T^f'  ,  'T^n');
plot(0:dt:T,Tn,'b','linewidth',2)
lgnd=legend('T^n' );
set(lgnd,'FontSize',14,'Location', 'Best');
set(gca,'FontSize',16)
xlabel('time','FontSize',16) 
ylabel('info flow','FontSize',16) 

subplot(3,2,[5,6]);
plot(0:dt:T,mean(x),'r',...
    0:dt:T,mean(y),'b','linewidth',2)
lgnd=legend('\mu_x'  ,  '\mu_y');
set(lgnd,'FontSize',10,'Location', 'Best');
set(gca,'FontSize',16)
xlabel('time','FontSize',16) 
ylabel('mean ','FontSize',16) 


subplot(3,2,3);
plot(0:dt:T,R11,'r',...
    0:dt:T,R22,'blue','linewidth',2)
lgnd=legend('R11' ,'R22');
set(lgnd,'FontSize',10,'Location', 'Best');
set(gca,'FontSize',16)
xlabel('time','FontSize',16) 
ylabel('variance ','FontSize',16)

subplot(3,2,4);
plot(0:dt:T,R12,'g','linewidth',2)
lgnd=legend('R12');
set(lgnd,'FontSize',10,'Location', 'Best');
set(gca,'FontSize',16)
xlabel('time','FontSize',16) 
ylabel('covariance ','FontSize',16)

% subplot(3,2,3);
% plot(0:dt:T,R11,'r',...
%     0:dt:T,R12,'--b',...
%     0:dt:T,R22,'--g','linewidth',2)
% lgnd=legend('R11'  ,  'R12','R22');
% set(lgnd,'FontSize',14,'Location', 'Best');
% set(gca,'FontSize',16)
% xlabel('time','FontSize',16) 
% ylabel('covariance ','FontSize',16)

 figure();
time=[1,2,4,6];
for i=1:4
time1=time(i);
%
subplot(2,2,i);
scatter(x(:,time1/dt),y(:,time1/dt))
set(gca,'FontSize',16)
xlabel('x','FontSize',16) 
ylabel('y ','FontSize',16)
title([ 'T= ', num2str(time1), ',  \sigma_{22}= ', num2str(sigma22)],'FontSize',16) 
end