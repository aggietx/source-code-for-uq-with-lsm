function ene2d=plotest2d(Kmax,kk,eneest)
ene2d=zeros(2*Kmax+1,2*Kmax+1);

for i=-Kmax:Kmax
    for j=-Kmax:Kmax
      i1=i+Kmax+1;
      j1=j+Kmax+1;
      id=getind(i,j,kk);
      ene2d(j1,i1)=eneest(id);
      
    end
end