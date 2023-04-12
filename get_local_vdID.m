function idmeshdrift=get_local_vdID(nx,gamma,x1,y1,L1,indmodes)
n=nx;%%wavenumber
lcase=2;
if lcase==1
[xa,ya]=meshgrid(0:L1/(n-1):L1,0:L1/(n-1):L1);
else
kmax=(n-1)/2;
% % [xa,ya]=meshgrid(0:L1*2*kmax/n/(n-1):L1*2*kmax/n,0:L1*2*kmax/n/(n-1):L1*2*kmax/n);
[xa,ya]=meshgrid(-L1*kmax/n:L1*2*kmax/n/(n-1):L1*kmax/n,-L1*kmax/n:L1*2*kmax/n/(n-1):L1*kmax/n);
x1=x1-L1/2;
y1=y1-L1/2;
end

xx=xa(:);
yy=ya(:);
Lmesh=size(xx,1);

L=size(x1,1);
dx1 = abs(xx * ones(1,L) - ones(Lmesh,1) * x1'); 
dy1 = abs(yy * ones(1,L) - ones(Lmesh,1) * y1'); 
% dy=dy1;dx=dx1;

dx2 = abs(xx * ones(1,L) - ones(Lmesh,1) * x1' + L1); 
dx3 = abs(xx * ones(1,L) - ones(Lmesh,1) * x1' - L1); 
dy2 = abs(yy * ones(1,L) - ones(Lmesh,1) * y1' + L1); 
dy3 = abs(yy * ones(1,L) - ones(Lmesh,1) * y1' - L1);
% 
distance_x_temp = min(dx1, dx2); % the distance in the x direction
dx = min(distance_x_temp, dx3); % the distance in the x direction
distance_y_temp = min(dy1, dy2); % the distance in the y direction
dy = min(distance_y_temp, dy3); % the distance in the y direction



d=sqrt(dx.^2+dy.^2)-gamma;
idmeshdrift=cell(Lmesh,1);
% parfor i=1:Lmesh
ind=1;
for i=1:Lmesh    
%     if ismember(i,indmodes)
    tempid=find(d(i,:)<0);
%     if size(tempid,2)==0
%         disp('no local obs v')
%     end  
%     [xx(i),yy(i)]
%     [x1(tempid),y1(tempid)]
    idmeshdrift{i}=tempid;
%     pause()
ind=ind+1;
%     end
end
% idmeshdrift=idmeshdrift(indmodes);