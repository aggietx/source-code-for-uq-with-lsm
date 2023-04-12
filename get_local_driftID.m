function [driftid,driftidpos]=get_local_driftID(gamma,x1,y1,x_obs,y_obs)
%%% gamma is the raidus
L=size(x1,1);
x1=real(x1);
y1=real(y1);
dx1 = abs(x1 * ones(1,L) - ones(L,1) * x_obs'); 
dy1 = abs(y1 * ones(1,L) - ones(L,1) * y_obs'); 


dx2 = abs(x1 * ones(1,L) - ones(L,1) * x_obs' + 2*pi); 
dx3 = abs(x1 * ones(1,L) - ones(L,1) * x_obs' - 2*pi); 
dy2 = abs(y1 * ones(1,L) - ones(L,1) * y_obs' + 2*pi); 
dy3 = abs(y1 * ones(1,L) - ones(L,1) * y_obs' - 2*pi);

distance_x_temp = min(dx1, dx2); % the distance in the x direction
dx = min(distance_x_temp, dx3); % the distance in the x direction
distance_y_temp = min(dy1, dy2); % the distance in the y direction
dy = min(distance_y_temp, dy3); % the distance in the y direction
% % dy=dy1;dx=dx1;


d=sqrt(dx.^2+dy.^2)-gamma;

driftid=cell(L,1);
driftidpos=cell(L,1);
% parfor i=1:L
for i=1:L
    tempid=find(d(:,i)<0);
    if size(tempid,1)==0
%         fprintf('no local obs drift %d\n',i)
    end
    pos=find(tempid==i);
    driftidpos{i}=pos;
    driftid{i}=tempid;
end
