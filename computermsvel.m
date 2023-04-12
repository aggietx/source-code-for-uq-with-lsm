function [rmsvel]=computermsvel(uex,uey,uestx,uesty,printindex)
%%% no normalization
uex=real(uex);
uey=real(uey);
uestx=real(uestx);
uesty=real(uesty);
rmsvel=rmse(uestx(:),uex(:),0)/2+rmse(uesty(:),uey(:),0)/2;


if printindex
    fprintf('rms ( time) %2.3f \n ',rmsvel);
end

