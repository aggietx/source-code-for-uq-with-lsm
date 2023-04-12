function [redind,redkk]=getreducedind(kk,Kmax,Kuse)
redind=[];
for ix=-Kmax:Kmax
    for iy=-Kmax:Kmax
        if sqrt(ix^2+iy^2)<=Kuse
        i1=getind(ix,iy,kk);redind=[redind;i1];
        end
    end
end
redkk=kk(:,redind);
    