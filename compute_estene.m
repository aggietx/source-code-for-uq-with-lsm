% Nit=5;
% allene12=allene;save allene12 allene12
% allene24=allene;save allene24 allene24
% load allene12;

if 1
%%%% compare the estimate energy and exact energy
ndata=floor(sqrt(2)*K_max);
allene=zeros(ndata,Nit+2);%%%% 1 intital,iteration, exact
enetraj=zeros(ndata,1);%%%% ene compute with trajectory of data
for i=1:ndata
    ind=getindrange(kk,i,K_max);
% %     enetraj(i)=mean(Rexact(ind));
    for j=1:Nit+2
        tempsig=sqrt(2)*allsigmak2(:,j);
        temp=alldkfit(:,j);
    allene(i,j)=mean( tempsig(ind).^2./temp(ind)/2);
    end
end
end
fprintf('L is %d\n',L);
return
% % allene12=allene;save allene12 allene12
load allene6;allene6(:,1:end-1)=-allene6(:,1:end-1);
load allene12;allene12(:,1:end-1)=-allene12(:,1:end-1);
load allene24;allene24(:,1:end-1)=-allene24(:,1:end-1);
load allene48;allene48(:,1:end-1)=-allene48(:,1:end-1);
load allene72;allene72(:,1:end-1)=-allene72(:,1:end-1);
load allene96;allene96(:,1:end-1)=-allene96(:,1:end-1);
      figure()
plot(allene12(:,end),'*-b','linewidth',2);hold on;
plot(allene6(:,1+Nit),'*-r','linewidth',2);hold on;
% plot(allene12(:,1+Nit),'*-k','linewidth',2);hold on;
plot(allene24(:,1+Nit),'*-g','linewidth',2);hold on;
plot(allene48(:,1+Nit),'*-k','linewidth',2);hold on;
% plot(allene72(:,1+Nit),'*-g','linewidth',2);hold on;
plot(allene96(:,1+Nit),'*-c','linewidth',2);hold on;
% plot(edk(:,6),'*-k','linewidth',2);hold on;
lgnd=legend('Exact energy','L=6','L=24','L=48','L=96');
set(lgnd,'FontSize',16);
setgca(18)
xlabel('|k|','FontSize',18);
xlim([1,size(allene12,1)]);
% title('Relative L_2 error of d_k','FontSize',18);
% xlim([6,96])
% set(gca,'XTick',[6,12,24,48,72,96])
