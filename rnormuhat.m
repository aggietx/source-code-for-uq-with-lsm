function rnormuhat(gamma_mean_trace,u_hat,ix,iy,kk)
ind=getind(ix,iy,kk);
rnorm(real(gamma_mean_trace(ind,:)),real(u_hat(ind,:)));

