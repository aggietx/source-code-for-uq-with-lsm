% 
% nmax=floor(sqrt(2)*Kmax);
% enekmean=zeros(nmax,1);
% eneksum=zeros(nmax,1);
% for i=1:nmax
%     
%     ind=[];
% for ix=-Kmax:Kmax
%     for iy=-Kmax:Kmax
%        if sqrt(ix^2+iy^2)<=i && sqrt(ix^2+iy^2)>i-1
%          ind=[ind;getind(ix,iy,kk)];
%        end     
%     end
%    
% end
% enekmean(i)=sum(ene(ind))/length(ind);
% eneksum(i)=sum(ene(ind));
% end

load ene
% % ene=[enekmean,eneksum,enekmean1,eneksum1];%%%  simulate/analytic
% k_0=2;E_0=1;alpha=3;
kcoeff=k_0*E_0*k_0^alpha;
alpha=3;
figure()
k=(1:size(ene,1));
plot(k,(ene(:,3)),'*-r','linewidth',2)
% hold on;plot(k(3:end),(k(3:end).^(-alpha)),'*-b','linewidth',2)
 setgca(16)
 
a=ene(:,1)/ene(3,1);
b=ene(:,3)/ene(3,3);
figure();plot(a,'*-r');hold on;plot(b,'*-g')