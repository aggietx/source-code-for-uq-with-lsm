function omegak=formomega(dimuhat,omega,kk,Kmax)
omegak=omega*ones(dimuhat,1);
omegak(1)=omega;
for iy=-Kmax:Kmax;
    for ix=0:Kmax;
%         ix=3;iy=2;
        i1=getind(ix,iy,kk);i2=getind(-ix,-iy,kk);
        omegak(i2)=-omegak(i1);
    end
end