%bus data generation
%adds another variable
%change from orig folder to go with Jon's parameters for Z and X from 3_z

function [Y,X,Z,Xstate,Zstate,State,FVT, p1]=genbus4_works(alpha,N,T,xtran,xtranc,xbin,zbin,xval,zval);

Beta=alpha(4);
Pi=alpha(5);
eul=0.577215665;
FV=zeros(xbin*zbin,2,T+10+1);

FV=V(alpha,T,xtran,xbin,zbin,xval);

Y=zeros(N,T+10);
X=zeros(N,1);
Zstate=ceil(length(zval)*rand(N,1));
State=(rand(N,1)>Pi);
p1 = ones(N,T+1+10);
Xstate=ones(N,T+1+10);
Draw=rand(N,T+10);
Draw2=rand(N,T+10);
%Draw3=rand(N,T+10);
Xstate=X+1;
for n=1:N;
    
    Z(n,1)=zval(Zstate(n,1));
    adj0=(Zstate(n)-1)*xbin+1;
    z2=(Zstate(n)-1)*xbin+1;
    z3=z2+xbin-1;
    for t=1:T+10;
        adj=Xstate(n,t)+(Zstate(n)-1)*xbin;
        FVT(n,t)=(xtran(adj,:)*FV(z2:z3,State(n)+1,t+1))-(xtran(adj0,:)*FV(z2:z3,State(n)+1,t+1));
        util1=alpha(1)+alpha(2)*X(n,t)+alpha(3)*State(n) + (xtran(adj,:)*FV(z2:z3,State(n)+1,t+1));
        util0=(xtran(adj0,:)*FV(z2:z3,State(n)+1,t+1));
        dem=exp(util1)+exp(util0);
        p0=exp(util0)./dem;
        p1(n, t) = 1-p0;
        Y(n,t)=1-(Draw(n,t)<p0);
        Xstate(n,t+1)=1+(Y(n,t)==1).*sum((Draw2(n,t)*ones(1,xbin-1))>xtranc(Xstate(n,t),1:xbin-1,Zstate(n,1)),2);
        X(n,t+1)=xval(Xstate(n,t+1));
    end;
end;

%Y=Y(:,11:T+10);
%X=X(:,11:T+10);
%Xstate=Xstate(:,11:T+10);
%FVT=FVT(:,11:T+10);

        

    
