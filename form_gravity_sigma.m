function sigma=form_gravity_sigma(kk0,sigma0)
dimuhat0=length(kk0);
sigma=zeros(2*dimuhat0,2*dimuhat0);

for i=1:length(kk0)
    sigma_g=sigma0(i);
ix=kk0(1,i);iy=kk0(2,i);[a]=find(kk0(1,:)==-ix);[b]=find(kk0(2,:)==-iy);
j=intersect(a+dimuhat0,b+dimuhat0); 
        sigma(i,i) = 1 / sqrt(2) * sigma_g;
        sigma(j, j) = -1i / sqrt(2) * sigma_g;
        sigma(i, j) = 1i / sqrt(2) * sigma_g;
        sigma(j, i) = 1 / sqrt(2) * sigma_g;
end
% sigma(1,1)=sigma_0(1);
% sigma(1+dimuhat0,1+dimuhat0)=1i*sigma_g(1);