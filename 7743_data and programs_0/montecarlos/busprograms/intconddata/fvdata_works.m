%calculates fv terms given the data 

function [fvt1, FV1]=fvdata(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T);


    FV1=zeros(tbin,2,T+1);
    
    for t=2:T;
        for s=0:1;
            FV1(:,s+1,t)= log(1+exp(-kron([1 t/10 (t/10)*(t/10)],[RX1 s*RX1])*b1));
        end;
    end;

    for n=1:N;
        
        adj0=(Zstate(n)-1)*xbin+1;
        z2=(Zstate(n)-1)*xbin+1;
        z3=z2+xbin-1;
        
        for t=1:T;
            adj=z2+Xstate(n,t)-1; 
            for s=1:2;                           
                FVT1(n,t,s)=(xtran(adj,:)-xtran(adj0,:))*FV1(z2:z3,s,t+1);
            end;
        end;
        
    end;
    fvt1=reshape(FVT1,N*T*2,1);