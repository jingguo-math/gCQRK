function [U,E]=cqrk_varn0_genkTrap(sol,f,Kfun,RK,N,Tf,alpha,grad)
% Evaluate K(partial_t)f at prescribed general time mesh tn, n=1,...,N
% with K(z)=1/(z^\alpha+1)
%using quadrature on the hyperbola to compute the history part
%It uses gCQ based on a Runge--Kutta method
%ref1: Generalized convolution quadrature with variable time
% stepping. Part II: Algorithm and numerical results
%=============== Modified by Jing Guo, 15-05-2024===============
%%
tolquad=1e-14;% tolerance of the qudrature rule
[s,A,c,b,eiginvA,V,invV]=setsolver(RK); %Buther tableau and eigenvalue decomposition of A^(-1)
eigA=eig(A);
est=[zeros(1, s-1) 1];
I=speye(s);
invA=A\I;
Mat1=zeros(s);
Mat1(:,end)=ones(s,1);
if RK == 1
    R = @(z) 1./(1 - z);
elseif RK == 2
    R = @(z) (2*z + 6)./(z.^2 - 4*z + 6);    
else
    R = @(z) est * (V * diag(1 ./ (1 - z * eigA)) * invV) * ones(s, 1);
end
if grad>20
    n0=N;
else
    n0=min(5,N); %Direct steps for RK-based gCQ
end
tvec = ((0:N)*Tf^(1/grad)/N).^grad;
hvec = tvec(2:end)-tvec(1:end-1);
%Approximation of convolution at first n0 steps
U=zeros(s,N+1);
E=zeros(1,N);
%Parameters for elliptic quadrature
if RK==1
    circ_radi=2;%choice of the radius of the circle contour
else
    circ_radi=1;
end
qvec=zeros(1,N-n0);
if n0==1
    %First step Omega0 = Kfun(A^(-1)/dt)
    Omega0 = Kfun(invA/hvec(1));
    tn=hvec(1)*c';
    U(:,2) = Omega0*f(tn);
else %n0>1
    Omega0 = Kfun(invA/hvec(1));
    tn=hvec(1)*c;
    U(:,2) = Omega0*f(tn);
    %Parameters for elliptic quadrature
    lamda=max(abs(eiginvA));
    
    switch circ_radi
        case 1
            m=min(real(eiginvA))/max(hvec(1:n0));
            M=5*lamda/min(hvec(1:n0));
        case 2
            m=1/max(hvec(1:n0));
            M=1/min(hvec(1:n0));
    end
    qvec(1)=M/m;
    a1=max(1,log(M)/log(m)/2);
    nq=ceil(n0^a1*log(n0)*(log(n0)+log(1/tolquad)));%ref1,Page 10, Cor 16
    [Z,weights]=elliptic_quadrature(m,M,nq);
    %Valores de la transformada de Laplace en el contorno
    E(1)=U(end,2)-sol(tvec(2));
    Yz=zeros(s,length(Z));
    for ll=3:n0+1
        tl=tvec(ll-1)+c*hvec(ll-1);
        %integral along circular contour
        int=0;
        for mm=1:length(Z)
            Yz(:,mm)=(I-hvec(ll-2)*A*Z(mm))\(Yz(end,mm)*ones(s,1)+...
                hvec(ll-2)*A*f(tvec(ll-2)+c*hvec(ll-2)));
            int=int+weights(mm)*Kfun(Z(mm))*((I-hvec(ll-1)*Z(mm)*A)\(Yz(end,mm)*ones(s,1)));
        end
        U(:,ll)=real(Kfun(invA/hvec(ll-1))*f(tl)+int);
        E(ll-1)=U(end,ll)-sol(tvec(ll));
    end
end
%Hasta aqui parece estar bien.
%%
%Approximation of convolution for n>n0

%Quadrature for the memory term
if n0<N
    [X,W]=Trap_quadrature(alpha,400,400,3/40);
    X=X';
    W=W';
    %Correction factor for memory term
    factor=ones(1,length(X));
    RVn=zeros(s,s*length(X));
    rn=ones(1,length(X));
    for ll=1:n0
        for mm=1:length(X)
            factor(mm)=factor(mm)/R(-hvec(ll)*X(mm));
        end
    end
    %For ODEs compression of the memory
    Yx=zeros(s,length(X));
end

for n=n0+1:N
    
    %Computation of first n0 weights for Local Term
    hloc=hvec(n-n0+1:n);
    %Parameters for elliptic quadrature
    %         circ_radi=1;%choice of the radius of the circle contour
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
    Yz=zeros(s,length(Z));
    for ll=1:n0-1
        %integral along circular contour
        for mm=1:length(Z)
            Yz(:,mm)=(I-hloc(ll)*A*Z(mm))^(-1)*(Yz(end,mm)*ones(s,1)+hloc(ll)*A*f(tvec(n-n0+ll)+c*hloc(ll)));
        end
    end
    int=zeros(s,1);
    for mm=1:length(Z)
        int = int + weights(mm)*Kfun(Z(mm))*((I-hvec(n)*Z(mm)*A)\(Yz(end,mm)*ones(s,1)));
    end
    local=int;
    
    %Memory Term for n>n0
    for mm=1:length(X)
        factor(mm)=factor(mm).*R(-hvec(n-n0)*X(mm))/R(-hvec(n)*X(mm));%？？？？
    end
    %Update ODEs up to t_{n-n0} and quadrature
    memory=0;
    for mm=1:length(X)
        Yx(:,mm)=(I+hvec(n-n0)*X(mm)*A)\(ones(s,1)*Yx(end,mm)+hvec(n-n0)*A*f(tvec(n-n0)+c*hvec(n-n0)));
        RVn(:,(mm-1)*s+1:mm*s)=(I+hvec(n)*X(mm)*A)\Mat1;
        rn(mm)=R(-hvec(n)*X(mm));
        memory=memory+RVn(:,(mm-1)*s+1:mm*s)/rn(mm)*(W(mm)/factor(mm)*Yx(:,mm));
    end
    tn=tvec(n)+hvec(n)*c;
    U(:,n+1)=real(Kfun(invA/hvec(n))*f(tn)+local+memory);
    E(n)=U(end,n+1)-sol(tvec(n+1));
end
E=max(abs(E));
end
