p=0.5;Nit=3;
p10=p*10;
fname = sprintf('ccall_est1D_it%dp%d.bin', Nit,p10);
 fid=fopen(fname,'r');ccall=fread(fid,'single');fclose(fid);%%% read
 fname = sprintf('rmsall_est1D_it%dp%d.bin', Nit,p10);
 fid=fopen(fname,'r');rmsall=fread(fid,'single');fclose(fid);%%% read