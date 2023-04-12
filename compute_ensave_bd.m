function [avexens,xensnew1]=compute_ensave_bd(xensnew,bdryh,MM)
% MM=size(xensnew,2);
% L=size(xensnew,1);
% idinnerx=find(sum(abs(xensnew-pi)<pi-bdryh,2)==MM);
% % idinnery=find(sum(abs(yensnew-pi)<pi-bdryh,2)==MM);
% idbdx=setdiff(1:L,idinnerx);
idbdx=find(sum(abs(xensnew-pi)<pi-bdryh,2)<MM);

% leftid=find(sum(abs(xensnew-bdryh/2)<bdryh/2,2)>=sum(sum(abs(xensnew-2*pi)<bdryh/2,2)));
% rightid=find(sum(abs(xensnew-bdryh/2)<bdryh/2,2)>=sum(sum(abs(xensnew-2*pi)<bdryh/2,2)));
% idbdy=setdiff(1:L,idinnery);
xensnew1=xensnew;
% yensnew1=yensnew;
 xensnew1(idbdx,:)=(2*pi - xensnew(idbdx,:)> pi) .* xensnew(idbdx,:) +...
       (2*pi - xensnew(idbdx,:) <= pi) .* (xensnew(idbdx,:) - 2*pi); 
%   yensnew1(idbdy,:)=(2*pi - yensnew(idbdy,:)> pi) .* yensnew(idbdy,:) +...
%        (2*pi - yensnew(idbdy,:) <= pi) .* (yensnew(idbdy,:) - 2*pi); 
avexens=sum(xensnew1,2)/MM;
% aveyens=sum(yensnew1,2)/MM;  