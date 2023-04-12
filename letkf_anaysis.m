imeshid=j;ndrift=length(idmeshdrift{imeshid});
imesh=j;
    localindx=idmeshdrift{imeshid};
xobslocal{imeshid}=temp_obs(localindx);yobslocal{imeshid}=temp_obs(localindx+nfloe);
localaveVx{imesh}=avexens(localindx);localVx{imesh}=V(localindx,:);
localaveVy{imesh}=aveyens(localindx);localVy{imesh}=V(localindx+nfloe,:);

      invR=diag((1/sigma_obs^2)*ones(ndrift*2,1));
      Ylf=[localVx{imesh};localVy{imesh}];
    aveylf=[localaveVx{imesh};localaveVy{imesh}];
    yol=[xobslocal{imesh};yobslocal{imesh}];
%     Q=(MM-1)*speye(MM)/(1+r)+1/(sigma_obs.^2)*(Ylf'*Ylf);
    Q=(MM-1)*speye(MM)/(1+r)+(Ylf'*invR*Ylf);

%     Pla=pinv(Q);
        [V1, D] = eig(Q/2+Q'/2); 
%     [V1, D] = eig(Q); 
    Pla=V1*diag(1./diag(D))*V1';%%%a=Pla-pinv(Q);max(abs(a(:)))  
% %     Wlatemp=sqrt(MM-1)*V1*diag(sqrt(1./diag(D)));
    Wlatemp=sqrt(MM-1)*V1*diag(sqrt(1./diag(D)))*V1';
%     avewla=1/(sigma_obs.^2)*Pla*Ylf'*(yol-aveylf);
      avewla=Pla*Ylf'*invR*(yol-aveylf);
% % % s=V*diag(sqrt(1./diag(D)));a=s*s'-Pla;max(abs(a(:)))
    Wla=Wlatemp+avewla*ones(1,MM);  
    ensuhatup([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0],:)=aveuhat([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0])*ones(1,MM)+...
        uhatmesh([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0],:)*Wla;
%     aveuhatup(imesh)=aveuhat(imesh)+Uhat(imesh,:)*avewla;
%     pause()

 if length(idmeshdrift{j})==0%%%% no update
     ensuhatup([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0],:)= uhatmesh([imesh;imesh+1*dimuhat0;imesh+2*dimuhat0],:);
 end
