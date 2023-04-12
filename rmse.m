function e=rmse(a,b,printindex)
a=a(:);
b=b(:);
n=length(a);
e=sqrt( sum((a-b).^2)/n);
if printindex==1
fprintf('rmse is %2.4f\n',e);
end