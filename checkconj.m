function checkconj(data,K_max,kk)
for ix=0:K_max
    for iy=0:K_max
 [~,s1]=getind(ix,iy,kk);
 
 if max(abs(imag(sum(data(s1,:)))))>10^(-6)
     disp('no conjuate')
 end
 if max(abs(real(data(s1(1),:))-real(data(s1(2),:)))) >10^(-6)
     disp('no conjuate, real are not the same')
 end
    end
end