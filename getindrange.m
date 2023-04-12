function ind=getindrange(kk,k1,K_max)
ind=[];
for i=-K_max:K_max
    for j=-K_max:K_max

      if norm([i,j])<=k1 && norm([i,j])>k1-1
          ind=[ind,getind(i,j,kk)];
      end        
    end
end