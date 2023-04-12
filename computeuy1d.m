function uy1d=computeuy1d(psik,Kmax,x)
kk0=(1:Kmax)';
a=psik;b=conj(a);%real(a)-1i*imag(a);
a1=1i.*a.*kk0;b1=-1i*b.*kk0;
uy1d = exp(1i * x * kk0')*a1+exp(-1i * x * kk0')*b1;
