%%%% prediction with nn
%% data
clear;load data_nntest
ix=1;iy=1;ind=getind(ix,iy,kk);
dataexact=real(u_hat(ind,:));
dataapp=real(gamma_mean_trace(ind,:));
diff=dataexact-dataapp;
nd=2000;ts=1000;
input=dataapp(ts:ts+nd);
output=diff(ts+1:ts+nd+1);

inputtest=dataapp(ts+nd+100:ts+nd+100+nd);
outputtest=diff(ts+nd+100+1:ts+nd+100+nd+1);
fprintf('# of sample is %d\n',nd);
%% nn
nlayer=1;
fprintf('nlayer is %d\n',nlayer);
% net = feedforwardnet(nlayer);
net=newff(input,output,nlayer);
%% training
disp('traning....')
Epochs=2000;
net.performFcn='sse';
net.trainParam.goal=0.001;%训练sse目标
net.trainParam.min_grad=1e-20;%最小梯度
net.trainParam.show = 200;
net.trainParam.epochs = Epochs;
net.trainParam.mc = 0.95;%动量因子
net.trainParam.lr=0.01;%学习率
tic;[net,tr] = train(net,input,output);toc

%% prediction
netout=net(input);
rnorm(netout,output); 
traincomparedata=[netout;output];save traincomparedata traincomparedata
rnorm(net(inputtest),outputtest); 
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