function sigmauhat=form_sigmatrix1d(dimuhat,sigmak2,kk,Kmax)
%%% form conjuate sigma matrix
sigmauhat = zeros(dimuhat, dimuhat);

for iy=0
    for ix=1:Kmax
%         ix=3;iy=2;
        i1=getind(ix,iy,kk);i2=getind(-ix,-iy,kk);
  sigma_B  =sigmak2(i1);    
% sigma_B=sqrt(sigmak2(i1));
  
%         sigmauhat(i1,i1) = sigma_B/sqrt(2);
%         sigmauhat(i2, i2) = -sqrt(-1)  * sigma_B/sqrt(2);
%         sigmauhat(i1, i2) = sqrt(-1)  * sigma_B/sqrt(2);
%         sigmauhat(i2, i1) = sigma_B/sqrt(2);
        sigmauhat(i1,i1) = sigma_B;
        sigmauhat(i2, i2) = -sqrt(-1)  * sigma_B;
        sigmauhat(i1, i2) = sqrt(-1)  * sigma_B;
        sigmauhat(i2, i1) = sigma_B;           
    end
end
sigmauhat(1,1)=sigmak2(1);