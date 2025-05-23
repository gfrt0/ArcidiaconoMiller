%bus data generation
%adds another variable
%change from orig folder to go with Jon's parameters for Z and X from 3_z

function [Y,X,Z,Xstate,Zstate,State,FVT,Adj]=genbus(alpha,N,T,xtran,xtranc,xbin,zbin,xval,zval,T2);

Beta=alpha(4);
Pi=alpha(5);
eul=0.577215665;
FV=zeros(xbin*zbin,2,T+10+1);

Adj=zeros(T+11,1);

t=2;

while t<T+12
    
    Adj(t)=.7*Adj(t)+.5*norminv(rand());
    t=t+1;
end

t=T+10;
while t>1
    for s=0:1;
        for z=1:zbin;
            for x=1:xbin
              adj=x+(z-1)*xbin;  
              util1=alpha(1)+Adj(t)+alpha(2)*xval(x)+alpha(3)*s+ xtran(adj,:)*FV((z-1)*xbin+1:z*xbin,s+1,t+1);
              util0=xtran(1+(z-1)*xbin,:)*FV((z-1)*xbin+1:z*xbin,s+1,t+1);
              FV(adj,s+1,t)=Beta*(log(exp(util1)+exp(util0)) + eul);
              
            end;   
        end;
    end;
    t=t-1;
end; 

Y=zeros(N,T2+10);
X=zeros(N,1);
Zstate=ceil(length(zval)*rand(N,1));
State=(rand(N,1)>Pi);

Xstate=ones(N,T2+1+10);
Draw=rand(N,T2+10);
Draw2=rand(N,T2+10);

for n=1:N;
    
    Z(n,1)=zval(Zstate(n,1));
    for t=1:T2+10;
        FVT(n,t)=(xtran(Xstate(n,t)+(Zstate(n)-1)*xbin,:)*FV((Zstate(n)-1)*xbin+1:Zstate(n)*xbin,State(n)+1,t+1))-(xtran((Zstate(n)-1)*xbin+1,:)*FV((Zstate(n)-1)*xbin+1:Zstate(n)*xbin,State(n)+1,t+1));
        util1=alpha(1)+Adj(t)+alpha(2)*X(n,t)+alpha(3)*State(n) + (xtran(Xstate(n,t)+(Zstate(n)-1)*xbin,:)*FV((Zstate(n)-1)*xbin+1:Zstate(n)*xbin,State(n)+1,t+1));
        util0=(xtran((Zstate(n)-1)*xbin+1,:)*FV((Zstate(n)-1)*xbin+1:Zstate(n)*xbin,State(n)+1,t+1));
        dem=exp(util1)+exp(util0);
        p0=exp(util0)./dem;
        Y(n,t)=1-(Draw(n,t)<p0);
        Xstate(n,t+1)=1+(Y(n,t)==1).*sum((Draw2(n,t)*ones(1,xbin-1))>xtranc(Xstate(n,t),1:xbin-1,Zstate(n,1)),2);
        X(n,t+1)=xval(Xstate(n,t+1));
    end;
end;

Y=Y(:,11:T2+10);
X=X(:,11:T2+10);
Xstate=Xstate(:,11:T2+10);
FVT=FVT(:,11:T2+10);