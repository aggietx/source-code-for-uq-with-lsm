function blowuptest(data,tol)

if max(abs(real(data(:))))>tol
  disp('blow up');
end