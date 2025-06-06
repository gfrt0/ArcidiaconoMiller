%shell for bus program
rng(44042, 'twister');

cd('D:/Dropbox (Personal)/git/ArcidiaconoMiller/7743_data and programs_0/montecarlos/busprograms/intconddata')

clear all;

Bccp=[];
Lccp=[]; 
Tccp=[];
Bf1=[];
Lf1=[]; 
Tf1=[];
Iccp=[];
Binit=[];

%parameter values
alpha=[2;-.15;1;.9;.4];
Beta=alpha(4);
eul=0.577215665;

T=20;
N=1000;

%optimization choices
%derivative is used in the minimization only for reduced form ccp's

o1=optimset('LargeScale','off','Display','off');
o2=optimset('LargeScale','off','Display','off','GradObj','on');
tol=.0000001;

%Create transition matrices;

zval=(.25:.01:1.25)';
zbin=length(zval);
xval=(0:.125:25)';
xbin=length(xval);
for z=1:zbin;
    [xtran(1+(z-1)*xbin:z*xbin,:),xtranc(:,:,z)]=xgrid(zval(z),xval);
end;

tbin=xbin*zbin;

%z and x values for each state

zvalr=kron(zval,ones(xbin,1));
xvalr=kron(ones(zbin,1),xval)./10;

%data for reduced form logits

RX1=[ones(zbin*xbin,1) xvalr zvalr xvalr.*zvalr xvalr.*xvalr zvalr.*zvalr];

%monte carlos

%starting values for FIML and CCP
alphaf=[2.2577;-.1384;.4;2.656;0;-.1;0];
alphac=[2.233;-.1339;.4;.9115];
alphac2=[0;-.1;0];

%generating the data
[Y,X,Z,Xstate,Zstate,State,TFV]=genbus4(alpha,N,T,xtran,xtranc,xbin,zbin,xval,zval);

y2=reshape(Y,N*T,1);
y2=[y2;y2];
x2=reshape(X,N*T,1)./10;
x2=[x2;x2];
z2=kron(ones(T,1),Z);
z2=[z2;z2];
s2=[zeros(N*T,1);ones(N*T,1)];
t2=kron([1:T]',ones(N,1))./10;
t2=[t2;t2];
stemp=[zeros(N,1);ones(N,1)];

%estimating with data ccps
fprintf('CCP estimation \n');

%tic

%setting up data for reduced form logit

xx=[ones(2*N*T,1) x2 z2 x2.*z2 x2.*x2 z2.*z2 s2 s2.*x2 s2.*z2 s2.*x2.*z2 s2.*x2.*x2 s2.*z2.*z2];
xx=[xx ((t2)*ones(1,12)).*xx ((t2.*t2)*ones(1,12)).*xx];

PType=.5*ones(2*N*T,1);
oPType=zeros(2*N*T,1);
Pi2=[.5 .5];

%estimating reduced form logit

b1=zeros(size(xx,2),1);

[b1]=fminunc('wlogitd',b1,o2,y2==0,xx,PType);

b1 = [
   -2.3537;
    2.5431;
   -0.2474;
    1.2398;
   -0.7335;
    0.1037;
    0.1215;
   -0.1541;
   -0.2770;
    0.1542;
    0.0451;
    0.1441;
   -0.5307;
    2.4721;
    1.5771;
   -1.9514;
   -0.3172;
   -0.5938;
   -0.2748;
    0.3611;
    0.6239;
   -0.3573;
   -0.1090;
   -0.3234;
    0.3161;
   -1.4664;
   -0.8205;
    0.8878;
    0.3427;
    0.2601;
    0.1200;
   -0.1596;
   -0.2722;
    0.1576;
    0.0487;
    0.1410
];

wlogitd(b1, y2 == 0, xx, PType)

% calculating fv terms

[fvt1, psi]=fvdata(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T);

%starting the EM algorithm
    
j=0;

bccp=alphac(1:4);
       
xccp=[ones(N*T*2,1) x2*10 s2];
intcondX=[ones(N,1) X(1:N,1) Z(1:N,1)];
binit=zeros(3,1);

cond=0;
lp=[];
while cond==0

    pi = round( mean(reshape(PType, N * T, 2)), 2);
    fprintf('EM iteration %d. Pi: %d\n', j, pi(1))

    %updating PType
    %first getting the type-specific likelihoods
    
    oPType=PType;
    Like=likeCCP(bccp,y2, [xccp fvt1]);
    Like2=reshape(Like,N,T,2);
    base=squeeze(prod(Like2,2));
    
    %now getting the initial condition parameters 
    
    [binit,lpx]=fminunc('intcond',binit,o1,base,intcondX);
    lp=[lp;lpx];
    
    %and the PType's
    
    [PType]=intcondP(binit,base,intcondX);
    
    PType=kron(ones(T,1),PType);

    PType=reshape(PType,(N*T*2),1);
    
    %estimating reduced form logit
        
    [b1]=fminunc('wlogitd',b1,o2,y2==0,xx,PType);
    
    %calculating fv terms
    
    [fvt1]=fvdata(b1,RX1,tbin,xbin,Zstate,Xstate,xtran,N,T);
    
    
    [bccp]=fminunc('wlogit',[bccp],o1,y2,[xccp fvt1],PType);
    
    %CHECKING CONVERGENCE
    
        if j>26;
            junk=abs((lp(j)-lp(j-25))./lp(j))<tol;
            junk2=abs((lp(j-1)-lp(j-26))./lp(j-1))<tol;
      
            cond=junk2.*junk;
            
            if j>1000;
                cond=1;
            end
        end        
            
    j=j+1;
end

tccp=toc

Tccp=[Tccp;tccp];
Bccp=[Bccp;[bccp;Pi2(2)]'];
Iccp=[Iccp;j];
Binit=[Binit;binit'];

% 
% MC 
% 
% save busdata0110 Bccp Tccp Iccp Binit Bf1 Tf1
% 
% 
