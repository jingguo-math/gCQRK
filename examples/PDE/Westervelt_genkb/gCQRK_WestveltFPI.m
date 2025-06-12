function [U]=gCQRK_WestveltFPI(u0,f,Kfun,alp,RK,N,Tf,grad,J,h,cl,v0,kap,c2d,r)

%the variable Runge-Kutta method for the linear wave equations with a
% convolution term:c2d(1-2*kapu)u_{tt}-cl\Deltau-kc*k\ast\Delta(u_t)=2*kap*(u_t)^2+f(x,t);
%It uses gCQ based on a Runge--Kutta method in time and compact FD in space
% the nonlinear system  is solved by the fixed point iteration
%ref1: Generalized convolution quadrature with variable time
%ref2: (Model) NUMERICAL ANALYSIS OF A TIME-STEPPING METHOD FOR
% THE WESTERVELT EQUATION WITH  TIME-FRACTIONAL DAMPING
%=============== Modified by Jing Guo, May, 2024===============
%% Initialization
tolquad=1e-14;% tolerance of the qudrature rule
[s,A,c,b,eiginvA,V,invV]=setsolver(RK); %Buther tableau and eigenvalue decomposition of A^(-1)
I=speye(s);
es=I(:,end);
invA=A\I;
rshp=@(x) reshape(x,s*(J-1),size(x,3));
rep=@(x) repmat(x,1,s);%replicate the matrix into (J-1)*s
tvec = ((0:N)*Tf^(1/grad)/N).^grad;
hvec = tvec(2:end)-tvec(1:end-1);
Lap=sparse(-2*diag(ones(1,J-1))+diag(ones(1,J-2),-1)+diag(ones(1,J-2),1))/(h^2);
Ms=sparse(10*diag(ones(1,J-1))+diag(ones(1,J-2),-1)+diag(ones(1,J-2),1))/12;
Lap=Ms\Lap; %fourth order finite differnce matrix
Ih=speye(J-1);
%%
if RK == 1
    R = @(z) 1./(1 - z);
elseif RK == 2
    R = @(z) (2*z + 6)./(z.^2 - 4*z + 6);
    qn= @(z) [9, 3-2*z]./(2*(z.^2-4*z+6));
    
else
    R = @(z) 1 + z*b'*((eye(s) - z*A)\ones(s,1)); % stability function
end
if grad>40
    n0=N;
else
    n0=min(5,N); %Direct steps for RK-based gCQ
end
%Parameters for elliptic quadrature
if RK==1
    circ_radi=2;%choice of the radius of the circle contour
else
    circ_radi=1;
end
%     [V D]=eig(A);
%     lamda=max(abs(diag(D)));
lamda=max(abs(eiginvA));
switch circ_radi
    case 1
        m=min(real(eiginvA))/max(hvec(1:n0));
        M=5*lamda/min(hvec(1:n0));
    case 2
        m=1/max(hvec(1:n0));
        M=1/min(hvec(1:n0));
end
a1=max(1,log(M)/log(m)/2);
nq=ceil(n0^a1*log(n0)*(log(n0)+log(1/tolquad)));%ref1,Page 10, Cor 16
[Z,weights]=elliptic_quadrature(m,M,nq);
dt=min(hvec(1:end-n0));%minimum time step in the history part
hmax=max(hvec(1:end-n0));%maximum time step in the history part
B=2;
Mat1=zeros(s);
Mat1(:,end)=ones(s,1);
%% compute numerical solution at t2
tolIt=10^-8;%tolerance of the fixed point iteration
if kap~=0
    ItFPI=500;
else
    ItFPI=0;
end
t1=tvec(1)+c'*hvec(1);

U=zeros(J-1,s,N+1);
U(:,:,1)=zeros(J-1,s);
Uin=rep(u0);%initial condition
L0=c2d*kron((invA/hvec(1))^2,Ih)-kron(I,cl*Lap)-kron(Kfun(invA/hvec(1))*(invA/hvec(1)),Lap);%coeffient matrix
rhs1=f(t1)+(c2d*rep(U(:,end,1))*transpose((invA/hvec(1))^2)+...
    c2d*rep(v0)*transpose(invA/hvec(1)))-(rep(Lap*U(:,end,1)))*transpose(Kfun(invA/hvec(1))*(invA/hvec(1)))+...
    cl*Lap*Uin+(Lap*rep(v0))*transpose(Kfun(invA/hvec(1)));%right hand side vector including the linear part
U_old=rshp(U(:,:,1));%initial condition of the fixed point iteration
%% fixed point iteration
itF=0;
while itF<=ItFPI
    U2=reshape(U_old,J-1,s);
    der2=(U2-rep(U(:,end,1)))*transpose((invA/hvec(1))^2)-rep(v0)*transpose(invA/hvec(1));
    rhsl=-2*c2d*kap*((U_old+rshp(Uin)).*rshp(der2));
    rhsr=2*kap*((U2-rep(U(:,end,1)))*transpose((invA/hvec(1)))).^2;
    rhs=rshp(rhs1+rhsr)-rhsl;
    U_new=L0\rhs;
%     norm(U_old-U_new,'fro')*h^0.5;
if max(max(abs(U_old-U_new)))<tolIt
      break
    end
    U_old=real(U_new);
    itF=itF+1;
end

U(:,:,2)=reshape(U_new,J-1,s);

%% compute numerical solution at tn 2<n<n0+1
Yz=zeros(J-1,s,length(Z));%solution of the ODE with z in the complex plane
for ll=2:n0
    int=0;
    for mm=1:length(Z)
        Yz(:,:,mm)=transpose((I-hvec(ll-1)*A*Z(mm))\transpose(((repmat(Yz(:,end,mm),1,s)+(U(:,:,ll)-repmat(U(:,end,ll-1),1,s))+...
            rep(v0)*transpose(hvec(ll-1)*A)))));
        int=int+weights(mm)*Kfun(Z(mm))*transpose((I-hvec(ll)*Z(mm)*A)\transpose(repmat(Yz(:,end,mm),1,s)));
    end
    int=real(int);
    tl=tvec(ll)+c'*hvec(ll);
    CV=U(:,:,ll-1:ll);
    CV=reshape(CV,s*(J-1),2);
    Ln=c2d*kron((invA/hvec(ll))^2,Ih)-kron(I,cl*Lap)-kron(Kfun(invA/hvec(ll))*(invA/hvec(ll)),Lap);%coefficent matrix
    rhs1= f(tl)+(c2d*rep(U(:,end,ll))*transpose((invA/hvec(ll))^2)+...
       c2d* ((varrk_1stder(hvec(ll-1),A,s,CV,J-1)*es)*ones(1,s))*transpose(invA/hvec(ll)))-...
        (rep(Lap*U(:,end,ll))*transpose(Kfun(invA/hvec(ll))*(invA/hvec(ll))))+...
        cl*Lap*Uin+ Lap*(int+rep(v0)*transpose(Kfun(invA/hvec(ll))));
    U_old=rshp(U(:,:,ll));
%% fixed point iteration
    itF=0;
    while itF<=ItFPI
        U2=reshape(U_old,J-1,s);
        der2=(U2-rep(U(:,end,ll)))*transpose((invA/hvec(ll))^2)-...
            ((varrk_1stder(hvec(ll-1),A,s,CV,J-1)*es)*ones(1,s))*transpose(invA/hvec(ll));
        rhsl=-2*c2d*kap*((U_old+rshp(Uin)).*rshp(der2));
        rhsr=2*kap*((U2-rep(U(:,end,ll)))*transpose((invA/hvec(ll)))).^2;
        rhs=rshp(rhs1+rhsr)-rhsl;
        U_new=Ln\rhs;
        norm(U_old-U_new,'fro')*h^0.5;
%         if norm(U_old-U_new,'fro')*h^0.5<tolIt
if max(max(abs(U_old-U_new)))<tolIt
            break
        end
        U_old=real(U_new);
        itF=itF+1;
    end
    U(:,:,ll+1)=reshape(U_new,J-1,s);
   
end

if n0<N
    [X,W,Q1,Q2s]=quadrature_cqw_grad(alp,tolquad,n0,tvec(end-n0),dt,hmax,B,RK);
    X=X';
    X=X+r;
    W=W';
    %Correction factor for memory term
    factor=ones(1,length(X));
    RVn=zeros(s,s*length(X));
    rn=ones(1,length(X));
    for ll=1:n0
        factor=factor./R(-hvec(ll)*X);
        
    end
    %For ODEs compression of the memory
    Yx=zeros(J-1,s,length(X));
end
%% compute numerical solution at tn n0+1<n<N
for n=n0+1:N
    %Computation of first n0 weights for Local Term
    hloc=hvec(n-n0+1:n);
%Parameters for elliptic quadrature
    switch circ_radi
        case 1
            m=min(real(eiginvA))/max(hloc);
            M=5*lamda/min(hloc);
        case 2
            m=1/max(hloc);
            M=1/min(hloc);
    end
    qvec(n-n0+1)=M/m;
    a1=max(1,log(M)/log(m)/2);
    nq=ceil(n0^a1*log(n0)*(log(n0)+log(1/tolquad)));%ref1,Page 10, Cor 16
    [Z,weights]=elliptic_quadrature(m,M,nq);
    %Valores de la transformada de Laplace en el contorno
    %     K = Kfun(Z);
    Yz=zeros(J-1,s,length(Z));
    for ll=1:n0-1
        %integral along circular contour
        int=0;
        for mm=1:length(Z)
            Yz(:,:,mm)=transpose((I-hloc(ll)*A*Z(mm))\transpose(repmat(Yz(:,end,mm),1,s)+U(:,:,n-n0+ll+1)-repmat(U(:,end,n-n0+ll),1,s)+...
                rep(v0)*transpose(hloc(ll)*A)));
            if ll==n0-1
                int = int + weights(mm)*Kfun(Z(mm))*transpose((I-hvec(n)*Z(mm)*A)\transpose(repmat(Yz(:,end,mm),1,s)));
            end
        end
    end
    
    
    local=int;
    memory=0;
    factor=factor.*(R(-hvec(n-n0)*X)./R(-hvec(n)*X));
    %Update ODEs up to t_{n-n0} and quadrature
    
    for mm=1:length(X)
        Yx(:,:,mm)=transpose((I+hvec(n-n0)*X(mm)*A)\transpose(repmat(Yx(:,end,mm),1,s)+U(:,:,n-n0+1)-repmat(U(:,end,n-n0),1,s)+...
            rep(v0)*transpose(hvec(n-n0)*A)));
        RVn(:,(mm-1)*s+1:mm*s)=(I+hvec(n)*X(mm)*A)\Mat1;
        rn(mm)=R(-hvec(n)*X(mm));
        memory=memory+transpose((RVn(:,(mm-1)*s+1:mm*s)/rn(mm))*transpose((W(mm)/factor(mm)*Yx(:,:,mm))));
    end
    tn=tvec(n)+c'*hvec(n);
    CV=U(:,:,n-1:n);
    CV=reshape(CV,s*(J-1),2);
    Ln=c2d*kron((invA/hvec(n))^2,Ih)-kron(I,cl*Lap)-kron(Kfun(invA/hvec(n))*(invA/hvec(n)),Lap);
    rhs1= f(tn)+(c2d*rep(U(:,end,n))*transpose((invA/hvec(n))^2)+...
        c2d*((varrk_1stder(hvec(n-1),A,s,CV,J-1)*es)*ones(1,s))*transpose(invA/hvec(n)))+...
        cl*Lap*Uin+ (Lap*rep(v0))*transpose(Kfun(invA/hvec(n)));
    rhs2=-(rep(Lap*U(:,end,n))*transpose(Kfun(invA/hvec(n))*(invA/hvec(n))))+Lap*(local+memory);
    rhs1=rhs1+rhs2;
    U_old=real(rshp(U(:,:,n)));
    %% fixed point iteration
    itF=0;
    while itF<=ItFPI
        U2=reshape(U_old,J-1,s);
        der2=(U2-rep(U(:,end,n)))*transpose((invA/hvec(n))^2)-...
            ((varrk_1stder(hvec(n-1),A,s,CV,J-1)*es)*ones(1,s))*transpose(invA/hvec(n));
        rhsl=-2*c2d*kap*((U_old+rshp(Uin)).*rshp(der2));
        rhsr=2*kap*((U2-rep(U(:,end,n)))*transpose((invA/hvec(n)))).^2;
        rhs=rshp(rhs1+rhsr)-rhsl;
        U_new=Ln\rhs;
        norm(U_old-U_new,'fro')*h^0.5;
        %     if norm(U_old-U_new,'fro')*h^0.5<tolIt
        if max(max(abs(U_old-U_new)))<tolIt
            break
        end
        U_old=real(U_new);
        itF=itF+1;
    end
    U(:,:,n+1)=reshape(U_new,J-1,s);
end
U=real(U+Uin);
end
