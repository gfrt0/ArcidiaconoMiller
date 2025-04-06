function FV=V(alpha,T,xtran,xbin,zbin,xval);

Beta=alpha(4);
eul=0.577215665;
FV=zeros(xbin*zbin,2,T+10+1);

t=T+10;
while t>1
    for s=0:1;
        for z=1:zbin;
            for x=1:xbin
              adj=x+(z-1)*xbin;  
              util1=alpha(1)+alpha(2)*xval(x)+alpha(3)*s                    + xtran(adj,:)*FV((z-1)*xbin+1:z*xbin,s+1,t+1);
              util0=                                                          xtran(1+(z-1)*xbin,:)*FV((z-1)*xbin+1:z*xbin,s+1,t+1);
              FV(adj,s+1,t)=Beta*(log(exp(util1)+exp(util0)) + eul);
              
            end;   
        end;
    end;
    t=t-1;
end; 