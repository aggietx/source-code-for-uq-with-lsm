function b=inan(a)
b= sum(isnan(a(:)));
if b>0
    disp('include nan')
end
c=sum(isinf(a(:)));
if c>0
    disp('include inf')
end