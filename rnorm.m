function e=rnorm(a,b)
e=norm(a(:)-b(:))/norm(b(:));
fprintf('relative error is %2.6f\n',e);
