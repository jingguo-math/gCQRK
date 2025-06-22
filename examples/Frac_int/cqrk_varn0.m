function [C,eamax]=cqrk_varn0(sol,f,Kfun,RK,N,Tf,alpha,grad)
% Evaluate I^{\alpha}[f] at prescribed general time mesh tn, n=1,...,N
%It uses gCQ based on a Runge--Kutta method
%ref1: Generalized convolution quadrature with variable time
% stepping. Part II: Algorithm and numerical results
%=====Written by Maria Lopez-Fernandez, Modified by Jing Guo, 23-04-2023====
%%
tolquad=1e-14;% tolerance of the qudrature rule
[s,A,c,b,eiginvA,V,invV]=setsolver(RK); %Buther tableau and eigenvalue decomposition of A^(-1)
% b=b';
% c=c';
I=eye(s); %estimated order of convergence plus one
if RK==1
    R=@(z) 1/(1-z);
else
    R=@(z) 1+z*b'*((I-z*A)\ones(s,1));%stability function
end

if grad>20
    n0=N;
else
    n0=min(5,N); %local steps for RK-based gCQ
end
B = 2;
tvec = ((0:N)*Tf^(1/grad)/N).^grad;
hvec = tvec(2:end)-tvec(1:end-1);
dt=min(hvec(1:end-n0));%minimum time step in the history part
hmax=max(hvec(1:end-n0));%maximum time step in the history part
%Approximation of convolution at first n0 steps
C=zeros(s,N);
%Parameters for elliptic quadrature
if RK==1
    circ_radi=2;%choice of the radius of the circle contour
else
    circ_radi=1;
end
if n0==1
    %First step Omega0 = Kfun(A^(-1)/dt)
    Omega0 = Kfun(eiginvA/hvec(1));
    Omega0 = V*diag(Omega0)*invV;
    tn=hvec(1)*c';
    C(:,1) = Omega0*f(tn);
else %n0>1
    Omega0 = Kfun(eiginvA/hvec(1));
    Omega0 = V*diag(Omega0)*invV;
    tn=hvec(1)*c;
    C(:,1) = Omega0*f(tn);
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
    
    a1=max(1,log(M)/log(m)/2);
    nq=ceil(n0^a1*log(n0)*(log(n0)+log(1/tolquad)));%ref1,Page 10, Cor 16
    [Z,weights]=elliptic_quadrature(m,M,nq);
    %Valores de la transformada de Laplace en el contorno
    K = Kfun(Z);
    U=zeros(s,length(Z));
    for ll=1:n0
        %integral along circular contour
        int=0;
        for mm=1:length(Z)
            U(end,mm)=R(hvec(ll)*Z(mm))*U(end,mm)+...
                hvec(ll)*b'*((eye(s)-hvec(ll)*Z(mm)*A)\f(tvec(ll)+c*hvec(ll)));
            int=int+weights(mm)*K(mm)*U(end,mm);
        end
        C(end,ll)=int;
    end
end
%Hasta aqui parece estar bien.
%%
%Approximation of convolution for n>n0

%Quadrature for the memory term
if n0<N
    [X,W]=quadrature_cqw_grad(alpha,tolquad,n0,tvec(end-n0),dt,hmax,B,RK);
    X=X';
    W=W';
    %Correction factor for memory term
    factor=ones(1,length(X));
    for ll=1:n0
        for mm=1:length(X)
            factor(mm)=factor(mm)/R(-hvec(ll)*X(mm));
        end
    end
    %For ODEs compression of the memory
    Y=zeros(s,length(X));
end

for n=n0+1:N
    %Computation of first n0 weights for Local Term
    if n0==1
        Omega0 = Kfun(eiginvA/hvec(n));
        Omega0 = V*diag(Omega0)*invV;
        local = Omega0*f(tvec(n+1));
    else %n0>1
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
        a1=max(1,log(m)/log(M)/2);
        nq=ceil(n0^a1*log(n0)*(log(n0)+log(1/tolquad)));
        [Z,weights]=elliptic_quadrature(m,M,nq);
        %Valores de la transformada de Laplace en el contorno
        K = Kfun(Z);
        U=zeros(s,length(Z));
        for ll=1:n0
            %integral along circular contour
            for mm=1:length(Z)
                U(end,mm)=R(hloc(ll)*Z(mm))*U(end,mm)+...
                    hloc(ll)*b'*((eye(s)-hloc(ll)*Z(mm)*A)^(-1)*f(tvec(n-n0+ll)+c*hloc(ll)));
            end
            local=0;
            for mm=1:length(Z)
                local=local+weights(mm).*K(mm).*U(end,mm);
            end            
        end
    end
    %Memory Term for n>n0
    for mm=1:length(X)
        factor(mm)=factor(mm).*R(-hvec(n-n0)*X(mm))/R(-hvec(n)*X(mm));%？？？？
    end
    %Update ODEs up to t_{n-n0} and quadrature
    memory=0;
    for kk=1:length(X)
        Y(end,kk)=R(-hvec(n-n0)*X(kk))*Y(end,kk)+...
            hvec(n-n0)*b'*((I+hvec(n-n0)*X(kk)*A)\f(tvec(n-n0)+c*hvec(n-n0)));
        memory=memory+W(kk)/factor(kk)*Y(end,kk);
    end
    C(end,n)=local+memory;
end
%Error
solvec = sol(tvec(2:end));
err = abs(C(end,:)-solvec);
eamax=max(err);
end
