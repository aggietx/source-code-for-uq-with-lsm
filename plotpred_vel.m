%%% compare exact and predict velocity
%%%% gauss and non gauss
close all;
mmplot=5;tlim=500;range=1:tlim;if testcase~=2;disp('error');end
figure();
step=500;
for j=1:3
    subplot(3,3,3*j-2)
uh=zeros(1+2*Kmax,1);
uh(1)=uens(1,step,j);
uh(2:2:end)=psikens(:,step,j);
uh(3:2:end)=conj(psikens(:,step,j));
 ux1=G1fixed*uh;uy1=G2fixed*uh;
if max(abs(imag( ux1(:))))>10^(-6) || max(abs(imag(uy1(:))))>10^(-6)
    disp('error, complex velicity')
else
    ux1=real(ux1);uy1=real(uy1);
end


% figure();
% subplot(2,1,1);
quiver(xx,yy,reshape(ux1,ny,nx),reshape(uy1,ny,nx),'k');
% hold on;scatter(x(:,end),y(:,end))
% title('exact velocity','fontsize',16);setgca(16)
if j==1
% title('truth','fontsize',16);   
title(['Truth at t= ', num2str(dt*step)],'fontsize',16)

end
xlim([0,L1]);ylim([0,L1]);setgca(16);
% pause(1)



    subplot(3,3,3*j-1)
uh=zeros(1+2*Kmax,1);
uh(1)=uens1(1,step,j);
uh(2:2:end)=psikens1(:,step,j);
uh(3:2:end)=conj(psikens1(:,step,j));
 ux1=G1fixed*uh;uy1=G2fixed*uh;
if max(abs(imag( ux1(:))))>10^(-6) || max(abs(imag(uy1(:))))>10^(-6)
    disp('error, complex velicity')
else
    ux1=real(ux1);uy1=real(uy1);
end


% figure();
% subplot(2,1,1);
quiver(xx,yy,reshape(ux1,ny,nx),reshape(uy1,ny,nx));
% hold on;scatter(x(:,end),y(:,end))
% title('exact velocity','fontsize',16);setgca(16)
if j==1
title('prediction (gauss)','fontsize',16);    
end
xlim([0,L1]);ylim([0,L1]);setgca(16);


    subplot(3,3,3*j)
uh=zeros(1+2*Kmax,1);
uh(1)=uens2(1,step,j);
uh(2:2:end)=psikens2(:,step,j);
uh(3:2:end)=conj(psikens2(:,step,j));
 ux1=G1fixed*uh;uy1=G2fixed*uh;
if max(abs(imag( ux1(:))))>10^(-6) || max(abs(imag(uy1(:))))>10^(-6)
    disp('error, complex velicity')
else
    ux1=real(ux1);uy1=real(uy1);
end


% figure();
% subplot(2,1,1);
quiver(xx,yy,reshape(ux1,ny,nx),reshape(uy1,ny,nx),'g');
% hold on;scatter(x(:,end),y(:,end))
% title('exact velocity','fontsize',16);setgca(16)
if j==1
title('prediction (non-gauss)','fontsize',16);    
end
xlim([0,L1]);ylim([0,L1]);setgca(16);


end
 filenameStrPdf = sprintf('velenscase%dtime%d.pdf',testcase,step);   
print(gcf, '-dpdf', '-r600', filenameStrPdf)

return
figure();
for id=1:16


subplot(4,4,id)
for iens=1:mmplot
    xtemp=fixtraj(xens(id,range,iens));
    ytemp=fixtraj(yens(id,range,iens));
plot(xtemp,ytemp,'k');hold on
end
plot(xtemp(1),ytemp(1),'or','linewidth',2);hold on
% xlim([0,2*pi]);ylim([0,2*pi])
%  xlim([-12,12]);ylim([-12,12])
% figure();

for iens=1:mmplot
    xtemp=fixtraj(xens2(id,range,iens));
    ytemp=fixtraj(yens2(id,range,iens));
plot(xtemp,ytemp,'b');hold on
end
% xlim([-12,12]);ylim([-12,12])
setgca(16)

for iens=1:mmplot
    xtemp=fixtraj(xens1(id,range,iens));
    ytemp=fixtraj(yens1(id,range,iens));
plot(xtemp,ytemp,'g');hold on
end
% xlim([-12,12]);ylim([-12,12])
setgca(16)
end
 filenameStrPdf = sprintf('trajenscase%dtlim%d.pdf',testcase,tlim);   
print(gcf, '-dpdf', '-r600', filenameStrPdf)
