halfindgb=[];halfindgb1=[];
for iy=-Kmax:Kmax;
    for ix=1:Kmax;
%         ix=3;iy=2;
        i1=getind(ix,iy,kk);
        halfindgb=[halfindgb;i1(1)];
                i1=getind(-ix,-iy,kk);
        halfindgb1=[halfindgb1;i1(1)];
    end
end

for iy=1:Kmax;
    for ix=0
%         ix=3;iy=2;
        i1=getind(ix,iy,kk);
        halfindgb=[halfindgb;i1(1)];
                        i1=getind(-ix,-iy,kk);
        halfindgb1=[halfindgb1;i1(1)];
    end
end

if length(halfindgb)~=(dimuhat0-1)/2
    disp('error, get half modes');
end

halfind=[halfindgb',2+dimuhat0:2*dimuhat0];halfind=halfind';
halfind1=[halfindgb1;negindkk(2:dimuhat0)+2*dimuhat0];
if norm(kk(:,halfind1)+ kk(:,halfind))~=0
  disp('error, get half modes seoncd part');  
end