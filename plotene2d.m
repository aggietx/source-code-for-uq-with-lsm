% load data_p1_L36_para2d
load data_p05_L36_para2d

ene2d=zeros(2*Kmax+1,2*Kmax+1);

for i=-Kmax:Kmax
    for j=-Kmax:Kmax
      i1=i+Kmax+1;
      j1=j+Kmax+1;
      id=getind(i,j,kk);
      ene2d(j1,i1)=eneest(id);
      
    end
end
figure();
imagesc(log(ene2d));colorbar;
setgca(18);
title('Recovered energy (log scale) with 2d model','fontsize',18)
axis off
return
load eneest1d
ene1d=ene2d(1+Kmax,1+Kmax:end);
figure();
subplot(3,1,1)
plot(1:Kmax+1,cumsum(eneexact)/sum(eneexact),'-*r','linewidth',2)
title('Exact energy','fontsize',18)
xlim([1,Kmax+1]);
setgca(18);
subplot(3,1,2)
plot(1:Kmax+1,cumsum(ene1d)/sum(ene1d),'-*r','linewidth',2)
title('Recovered energy with 2d model','fontsize',18)
xlim([1,Kmax+1]);
setgca(18);
% xlabel('|k|','fontsize',18);

subplot(3,1,3)
plot(1:Kmax+1,cumsum(eneest1d([1,2:2:end]))/sum(eneest1d([1,2:2:end])),'-*r','linewidth',2)
title('Recovered energy with 1d model','fontsize',18)
xlim([1,Kmax+1]);
setgca(18);
xlabel('|k|','fontsize',18);