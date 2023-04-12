function [rmsvel,ccvel]=rmscc(vapp,vex,printindex)
vapp=real(vapp);
vex=real(vex);
rmsvel=rmse(vapp(:),vex(:),0)/std(vex(:));
ccvel=corr2(vapp(:),vex(:));

if printindex
    fprintf('rms cc are %2.3f & %2.3f\n',rmsvel, ccvel);
end

