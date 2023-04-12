function  [enekmean,eneksum]=computeene1D(Kmax,ene,kk)
nmax=floor(sqrt(2)*Kmax);
enekmean=zeros(nmax,1);
eneksum=enekmean;
for i=1:nmax
    
    ind=[];
for ix=-Kmax:Kmax
    for iy=-Kmax:Kmax
       if sqrt(ix^2+iy^2)<=i && sqrt(ix^2+iy^2)>i-1
         ind=[ind;getind(ix,iy,kk)];
       end     
    end
   
end
enekmean(i)=sum(ene(ind))/length(ind);
eneksum(i)=sum(ene(ind));
end