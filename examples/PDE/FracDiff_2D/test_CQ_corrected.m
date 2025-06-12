%% Corrected BDF3 for 2d PDE D_t^\alpha u-c_u*Lap*u=f(f is a given function) at prescribed general time mesh tn, n=1,...,N 
 clear
%      RK=2;
Tf=1;%1(original)
alp=.75;
bet=alp;

J=256;



E=[];
c_u=1;
example=1;
lps=1; %left points
 a=-lps;
 b=lps;
 h=(b-a)/J;
 x1=((a+h):h:b-h)';
 y1=x1';
sol = @(t) cos(pi/2*x1)*cos(pi/2*y1).*t.^bet; 
f=@(t) cos(pi/2*x1)*cos(pi/2*y1)*gamma(bet+1).*t.^(bet-alp)/gamma(bet-alp+1)+...
     pi^2/2*cos(pi/2*x1)*cos(pi/2*y1).*c_u.*t.^bet;


%% 
% maxIt =4;
Nvec=2*2.^[1:7]'; 
maxIt=length(Nvec);
E= zeros(maxIt,1);
cput = zeros(maxIt,1);
%% Generate an initial mesh

%% 
for k = 1:maxIt
    
    tic
  [C,L,eL2max] =CQ_BDF3Corrected(sol,f,Nvec(k),Tf,J,alp,h,c_u);
    cput(k)=toc;
E(k,1)=eL2max;

end
E% Nvec
rate=log2(E(1:end-1)./E(2:end))
