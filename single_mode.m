% Set rng seed for reproducability.
rng(1)

% stepSize = 0.01; % Step size for numerical integration.
stepSize = 1/24; % Step size for numerical integration.
tFinal = 365;
tArray = 0:stepSize:tFinal;
nIterations = length(tArray);

% Setting true parameters.
dTruth = 3; % damping
phiTruth = 1; % imaginary component of damping
fTruth = 2; % forcing
sigTruth = 2; % stochastic forcing

% Initializing with NaN helps with debugging.
uArray = nan(1, nIterations);
uArray(1) = 0;

% Euler-Maruyama scheme.
for iIteration = 1:nIterations-1
    u = uArray(iIteration);
    
    % I prefer short variable names when implementing formulae.
    d = dTruth;
    phi = phiTruth;
    f = fTruth;
    sig = sigTruth;
    
    dt = stepSize;
    dW = sqrt(dt)*randn(1);
    
    du = ((-d + phi*1i)*u + f)*dt + sig*dW;
    
    uArray(iIteration+1) = u + du;
end

% uArray=uhaty(ind,:);
uMean = mean(uArray);
uVar = var(uArray);

% Typically, we would use autocorr from the econometrics toolbox. However
% it cannot handle complex variables. Instead we compute the
% autocorrelation using cross-correlation.
maxLag = 100;
% autoCorr = xcorr(uArray - uMean, maxLag, 'normalized');

autoCorr = xcorr(uArray - uMean);
autoCorr=autoCorr(length(uArray):end)/var(uArray,1)/length(uArray);

% xcorr calculates the correlation for positive and negative lags, but for
% the decorrelation time we only want positive lags. Then we integrate the
% autocorrelation function.
% uDecorr = trapz(stepSize, autoCorr(maxLag+1:end));
uDecorr = trapz(autoCorr(1:maxLag))*stepSize;

% Recovered parameters.
dRecovered = real(1/uDecorr);
phiRecovered = -imag(1/uDecorr);
fRecovered = uMean / uDecorr;
sigRecovered = sqrt(2*uVar/uDecorr);

% uArray=zeros(nIterations,1);
% uArray(1)=uhaty(ind,1);
% for iIteration = 1:nIterations-1
%     u = uArray(iIteration);
%     
%     % I prefer short variable names when implementing formulae.
%     d = dRecovered;
%     phi = phiRecovered;
%     f = fRecovered;
%     sig = sigRecovered;
%     
%     dt = stepSize;
%     dW = sqrt(dt)*randn(1);
%     
%     du = ((-d + phi*1i)*u + f)*dt + sig*dW;
%     
%     uArray(iIteration+1) = u + du;
% end
% plotpdf(real(uArray))
% plotpdf(real(uhaty(ind,:)))

fprintf('True damping: %f\n', dTruth);
fprintf('Recovered damping: %f\n', dRecovered);

fprintf('True phase: %f\n', phiTruth);
fprintf('Recovered phase: %f\n', phiRecovered);

fprintf('True forcing: %f\n', fTruth);
fprintf('Recovered forcing: %f\n', fRecovered);

fprintf('True noise: %f\n', sigTruth);
fprintf('Recovered noise: %f\n', sigRecovered);
