function x=fixtraj(x)
N=length(x);
tol=6;tol1=tol*2;tol2=tol*3;
for i=2:N
    x2=x(i);x1=x(i-1);
    if x2-x1>tol
        x(i)=x(i)-2*pi;
    elseif x2-x1<-tol
        x(i)=x(i)+2*pi;
    end
    
    x2=x(i);x1=x(i-1);
    if x2-x1>tol1
        x(i)=x(i)-2*pi*2;
    elseif x2-x1<-tol1
        x(i)=x(i)+2*pi*2;
    end   
%     
%         x2=x(i);x1=x(i-1);
%     if x2-x1>tol2
%         x(i)=x(i)-2*pi*3;
%     elseif x2-x1<-tol2
%         x(i)=x(i)+2*pi*3;
%     end  
%     x2=x(i);x1=x(i-1);
%     if x2-x1>6.2
%         x(i)=x(i)-2*pi;
%     end
%     if x2-x1<-6.2
%         x(i)=x(i)+2*pi;
%     end
    
end
return
if x(2)-x(1)>6
    x(1)=x(1)+2*pu;
end

if x(1)-x(2)<-6
    x(1)=x(1)-2*pu;
end
