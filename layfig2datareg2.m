clear;close all;
LL=[2,6,12,24,36];
rmsregII1=zeros(5,1);
for jj=1:5
   if jj==1 
load data_enkbf2_L3_p2
   elseif jj==2
load data_enkbf2_L6_p2
   elseif jj==3
load data_enkbf2_L12_p2
   elseif jj==4
load data_enkbf2_L24_p2
   elseif jj==5      
load data_enkbf2_L36_p2          
   end
estuhat1=estuhat;uhate1=uhate;estuhat1(2:end,:)=0;uhate1(2:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat1,uhate1,Kmax,K_max);
rmsregII1(jj,1)=meanrms1;
estuhat1=estuhat;uhate1=uhate;estuhat1(6:end,:)=0;uhate1(6:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat1,uhate1,Kmax,K_max);
rmsregII1(jj,2)=meanrms1;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat,uhate,Kmax,K_max);
rmsregII1(jj,3)=meanrms1;
end

rmsregII2=zeros(5,1);
for jj=1:5
   if jj==1 
load data_enkbf2_reduced_L3_p2 
   elseif jj==2
load data_enkbf2_reduced_L6_p2
   elseif jj==3
load data_enkbf2_reduced_L12_p2
   elseif jj==4
load data_enkbf2_reduced_L24_p2
   elseif jj==5      
load data_enkbf2_reduced_L36_p2        
   end
estuhat1=estuhat;uhate1=uhate;estuhat1(2:end,:)=0;uhate1(2:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat1,uhate1,Kmax,K_max);
rmsregII2(jj,1)=meanrms1;
estuhat1=estuhat;uhate1=uhate;estuhat1(6:end,:)=0;uhate1(6:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat1,uhate1,Kmax,K_max);
rmsregII2(jj,2)=meanrms1;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat,uhate,Kmax,K_max);
rmsregII2(jj,3)=meanrms1;
end

% load data_enkbf2_reduced_L12 
% 
rmsregII3=zeros(5,1);
for jj=1:5
   if jj==1 
load data1d_layer2_L3_p2
   elseif jj==2
load  data1d_layer2_L6_p2
   elseif jj==3
load data1d_layer2_L12_p2
   elseif jj==4
load  data1d_layer2_L24_p2
   elseif jj==5      
load data1d_layer2_L36_p2        
   end
estuhat1=estuhat;uhate1=uhate;estuhat1(2:end,:)=0;uhate1(2:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat1,uhate1,Kmax,K_max);
rmsregII3(jj,1)=meanrms1;
estuhat1=estuhat;uhate1=uhate;estuhat1(6:end,:)=0;uhate1(6:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat1,uhate1,Kmax,K_max);
rmsregII3(jj,2)=meanrms1;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat,uhate,Kmax,K_max);
rmsregII3(jj,3)=meanrms1;
end


rmsregII4=zeros(5,1);
for jj=1:5
   if jj==1 
load data1dreduced_layer2_L3_p2
   elseif jj==2
load  data1dreduced_layer2_L6_p2
   elseif jj==3
load data1dreduced_layer2_L12_p2
   elseif jj==4
load  data1dreduced_layer2_L24_p2
   elseif jj==5      
load data1dreduced_layer2_L36_p2       
   end
estuhat1=estuhat;uhate1=uhate;estuhat1(2:end,:)=0;uhate1(2:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat1,uhate1,Kmax,K_max);
rmsregII4(jj,1)=meanrms1;
estuhat1=estuhat;uhate1=uhate;estuhat1(6:end,:)=0;uhate1(6:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat1,uhate1,Kmax,K_max);
rmsregII4(jj,2)=meanrms1;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer1dfun(estuhat,uhate,Kmax,K_max);
rmsregII4(jj,3)=meanrms1;
end


rmsregII5=zeros(5,1);
for jj=1:5
   if jj==1 
load data2d_layer2_L3_p2
   elseif jj==2
load  data2d_layer2_L6_p2
   elseif jj==3
load data2d_layer2_L12_p2
   elseif jj==4
load  data2d_layer2_L24_p2
   elseif jj==5      
load data2d_layer2_L36_p2      
   end
estuhat1=estuhat;uhate1=uhate;estuhat1(2:end,:)=0;uhate1(2:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer2dfun(estuhat1,uhate1,Kmax);
rmsregII5(jj,1)=meanrms1;
ind=[1;getind(1,0,kk);getind(2,0,kk);getind(-1,0,kk);getind(-2,0,kk);];
ind1=setdiff(1:dimuhat,ind);
estuhat1=estuhat;uhate1=uhate;estuhat1(ind1,:)=0;uhate1(6:end,:)=0;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer2dfun(estuhat1,uhate1,Kmax);
rmsregII5(jj,2)=meanrms1;
[meanrmstime1,meanrms1,meancc1]=compute_errorlayer2dfun(estuhat,uhate,Kmax);
rmsregII5(jj,3)=meanrms1;
end
