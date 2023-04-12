function [dk,omegak,fk,sigmak,ene,sigmak0]=form_coeff(dimuhat,kk,E_0,k_0,alpha,Kmax,f0)
% f0=0;%%%% constant force
% E_0=1;
% k_0=2;
% alpha=3;
dk=zeros(dimuhat,1);
omegak=zeros(dimuhat,1);
fk=zeros(dimuhat,1);
sigmak0=zeros(dimuhat,1);
ene=zeros(dimuhat,1);
% d=1;nu=1;
d=0.3;nu=.05;
% d=0.0;nu=.05;
for i=2:dimuhat
   k=kk(:,i);ix=k(1);iy=k(2);
%         
%         ix=3;iy=2;i1=getind(ix,iy,kk);i2=getind(-ix,-iy,kk);s=[i1;i2];
% %         i3=getind(ix,-iy,kk);i4=getind(-ix,iy,kk);
% %         k=[ix,iy];
        if norm(k)<=k_0
            ene(i)=norm(k)*E_0;
        else
            ene(i)=k_0*E_0*abs(norm(k)/k_0)^(-alpha);
        end
       tempdk = d + nu*norm(k)^2;
       dk(i)=tempdk;
       fk(i)=f0;
       ek=ene(i);
       sigmak0(i) = sqrt(2*dk(i)*(2*ek-fk(i)^2/dk(i)^2));
end

% return
% % for iy=-Kmax:Kmax;
% %     for ix=0:Kmax;
% % %         ix=3;iy=2;
% %         i1=getind(ix,iy,kk);i2=getind(-ix,-iy,kk);
% %         omegak(i2)=-omegak(i1);
% %     end
% % end

sigmak = zeros(dimuhat, dimuhat);

for iy=-Kmax:Kmax;
    for ix=0:Kmax;
%         ix=3;iy=2;
        i1=getind(ix,iy,kk);i2=getind(-ix,-iy,kk);
         omegak(i2)=-omegak(i1);
         
  sigma_B  =sigmak0(i1);    
% sigma_B=sqrt(sigmak2(i1));
  
        sigmak(i1,i1) = sigma_B/sqrt(2);
        sigmak(i2, i2) = -sqrt(-1)  * sigma_B/sqrt(2);
        sigmak(i1, i2) = sqrt(-1)  * sigma_B/sqrt(2);
        sigmak(i2, i1) = sigma_B/sqrt(2);
           
    end
end
sigmak(1,1)=sigmak0(1);