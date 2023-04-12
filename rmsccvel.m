function [rmsvel,ccvel]=rmsccvel(uex,uey,uestx,uesty,printindex)
uex=real(uex);
uey=real(uey);
uestx=real(uestx);
uesty=real(uesty);
rmsvel=rmse(uestx(:),uex(:),0)/std(uex(:))/2+rmse(uesty(:),uey(:),0)/std(uey(:))/2;
ccvel=corr2(uestx(:),uex(:))/2+corr2(uesty(:),uey(:))/2;

if printindex
    fprintf('rms cc are %2.3f & %2.3f\n ',rmsvel, ccvel);
end

