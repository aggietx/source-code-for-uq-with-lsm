function ene=getene(data)
a=real(data);
b=imag(data);
ene=sum(a.^2+b.^2,2)/length(data)/2;