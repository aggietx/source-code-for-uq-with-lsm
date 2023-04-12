function [ind1]=getindpair(ix,iy,kk)
[a]=find(kk(1,:)==ix);
[b]=find(kk(2,:)==iy);
ind=intersect(a,b);

[a]=find(kk(1,:)==-ix);
[b]=find(kk(2,:)==-iy);
ind1=[ind,intersect(a,b)];