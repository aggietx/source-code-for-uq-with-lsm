%%%% prediction with lstm
%% data
clear;load data_nntest
ix=1;iy=1;ind=getind(ix,iy,kk);
dataexact=real(u_hat(ind,:));
dataapp=real(gamma_mean_trace(ind,:));
diff=dataexact-dataapp;
nd=200;ts=1000;
input=dataapp(ts:ts+nd);
output=diff(ts+1:ts+nd+1);

inputtest=dataapp(ts+nd+100:ts+nd+100+nd);
outputtest=diff(ts+nd+100+1:ts+nd+100+nd+1);
fprintf('# of sample is %d\n',nd);
%% nn
numFeatures=size(input,1);
numResponses=size(output,1);
numHiddenUnits=200;

layers = [
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    regressionLayer];
options=trainingOptions('adam',...
    'MaxEpochs',2000,...
    'GradientThreshold',1,...
    'InitialLearnRate',0.005,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod',125,...
    'LearnRateDropFactor',0.2,...
    'Verbose',0);
disp('training')
tic;net=trainNetwork(input,output,layers,options);toc

%% prediction
netout=predict(net, inputtest);
traincomparedata=[netout;output];save traincomparedata traincomparedata

rnorm(netout,output); 


return

x1=1:10;x2=2:11;y=3+x1*3+x2*4;x=[x1;x2];

net=newff(x,y,1);
net= train(net,x,y);
rnorm(net(x),y);



return
%%
gap2=gap;
ue=zeros(dimuhat,N);ue1=zeros(dimuhat,N);
ufree=zeros(dimuhat,N);
sigmafit=form_sigmatrix1(size(u_hat,1),allsigmak2(:,Nit),kk,Kmax);
disp('free run with fixed noise (exact parameter and estimated parameter) ..')
for i=2:N
     noise1=randn(dimuhat,1);
     noise2=randn(dimuhat,1);
    ue(:,i)=ue(:,i-1)+(dkexact+1i*omegaexact).*ue(:,i-1)*dt+...
        Fuexact*dt+sqrt(dt)*sigmauhat*noise1;
    ue1(:,i)=ue1(:,i-1)+(dkexact+1i*omegaexact).*ue1(:,i-1)*dt+...
        Fuexact*dt+sqrt(dt)*sigmauhat*noise2;
    ufree(:,i)=ufree(:,i-1)+(-alldkfit(:,Nit)+1i*allomegakfit(:,Nit)).*ufree(:,i-1)*dt+...
        allFfit(:,Nit)*dt+sqrt(dt)*sigmafit*noise2;
    if max(abs(real(ue(:,i))))>10^8 || max(abs(real(ufree(:,i))))>10^8
        disp('error, data blow up')
    end
%     tempcode=u_hat(:,i);
% %     utest1(:,N-i+1)=tempcode;
% %     utest(:,i)=tempcode;
    
end