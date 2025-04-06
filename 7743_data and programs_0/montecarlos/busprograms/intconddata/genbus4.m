%bus data generation
%adds another variable
%change from orig folder to go with Jon's parameters for Z and X from 3_z

function [Y,X,Z,Xstate,Zstate,State,FVT]=genbus(alpha,N,T,xtran,xtranc,xbin,zbin,xval,zval);

FV=V(alpha,T,xtran,xbin,zbin,xval);

Pi=alpha(5);
Xstate=ones(N, T+10);
Zstate=ceil(length(zval)*rand(N,1));
State=(rand(N,1)>Pi);
Y=zeros(N,T+10);
X=zeros(N,1);
Z=zeros(N,1);
p1 = ones(N,T+1+10);
Draw=rand(N,T+10);
Draw2=rand(N,T+10);
FVT = zeros(N, T+10);

for n=1:N;
    Z(n,1)=zval(Zstate(n,1));
    for t=1:T+10;
        Vslice = FV((Zstate(n)-1)*xbin+1:(Zstate(n))*xbin,State(n)+1,t+1);
        FVT(n,t)=(xtran(Xstate(n,t)+(Zstate(n)-1)*xbin,:) - xtran((Zstate(n)-1)*xbin+1,:)) * Vslice;
        deltav = alpha(1)+alpha(2)*X(n,t)+alpha(3)*State(n) + FVT(n,t);
        p0=1./(exp(deltav)+1);
        p1(n, t) = 1-p0;
        Y(n,t)=Draw(n,t)>p0;
        Xstate(n,t+1)=1+(Y(n,t)==1).*sum((Draw2(n,t)*ones(1,xbin-1))>xtranc(Xstate(n,t),1:xbin-1,Zstate(n,1)),2);
        X(n,t+1)=xval(Xstate(n,t+1));
    end;
end;

Y=Y(:,11:T+10);
X=X(:,11:T+10);
Xstate=Xstate(:,11:T+10);
FVT=FVT(:,11:T+10);