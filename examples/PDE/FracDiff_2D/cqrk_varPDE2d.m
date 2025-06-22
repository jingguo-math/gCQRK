function [C,eamax,eTf]=cqrk_varPDE2d(sol,f,RK,N,Tf,grad,J,alp,h,c_u)
% approximate the 2d PDE D_t^\alpha u-c_u*Lap*u=f(f is a given function) at prescribed general time mesh tn, n=1,...,N
%using quadrature on the hyperbola to compute the history part
%It uses gCQ based on a Runge--Kutta method in time and compact FD with
%discrete sine transform(refer to ref2) in space
%ref1: Generalized convolution quadrature with variable time
% stepping. Part II: Algorithm and numerical results
% ref2: Fast High-Order Compact Exponential Time Differencing Runge–Kutta Methods for Second-Order Semilinear Parabolic Equations
%=============== Modified by Jing Guo, May, 2024===============
%%


tolquad=1e-14;% tolerance of the qudrature rule
[s,A,c,b,eiginvA,V,invV]=setsolver(RK); %Buther tableau and eigenvalue decomposition of A^(-1)
I=speye(s);
Mat1=zeros(s);
Mat1(:,end)=ones(s,1);
if RK == 1
    R = @(z) 1./(1 - z);
elseif RK == 2
    R = @(z) (2*z + 6)./(z.^2 - 4*z + 6);
    qn= @(z) [9, 3-2*z]./(2*(z.^2-4*z+6));
    
else
  R = @(z) arrayfun(@(zi) 1 + zi * b' * ((eye(s) - zi * A) \ ones(s,1)), z); %stability function
end
if grad>20
    n0=N;
else
    n0=min(5,N); %Direct steps for RK-based gCQ
end
B = 2;
%graded for [0,Tf]
tvec = ((0:N)*Tf^(1/grad)/N).^grad;
hvec = tvec(2:end)-tvec(1:end-1);
dt=min(hvec(1:end-n0));%minimum time step in the history part
hmax=max(hvec(1:end-n0));%maximum time step in the history part
D=-4/h^2*(sin((1:(J-1))'*pi/(2*J))).^2;
D=c_u*(D+D'+1/12*(h^2+h^2)*D*D')./(1+1/12*(h^2*D+h^2*D'));%eigenvalue of the compact FD matrix,ref 2: Page 6
%Approximation of convolution at first N steps
C=zeros((J-1)^2,s,N+1);
C(:,:,1)=repmat(reshape(sol(0),(J-1)^2,1),1,s);
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
qvec=zeros(1,N-n0);
qvec(1)=M/m;
t2=c*hvec(1);
P1=((rf(f,t2,J,s))*transpose((hvec(1)*A)^alp))*transpose(invV);
P2=idstf(P1,J,s);
P3=1./(1-transpose(eig(A)).^alp.*(hvec(1)^alp*reshape(D,(J-1)^2,1)));
P3=P3.*P2;
C(:,:,2)=dstf(P3,J,s)*transpose(V);
E=[];
E=[E norm(real(reshape(C(:,end,2),J-1,J-1))-sol(tvec(2)),'fro')*h];
a1=max(1,log(M)/log(m)/2);
nq=ceil(n0^a1*log(n0)*(log(n0)+log(1/tolquad)));%ref1,Page 10, Cor 16
[Z,weights]=elliptic_quadrature(m,M,nq);
%Valores de la transformada de Laplace en el contorno

U=zeros((J-1)^2,s,length(Z));
for ll=3:n0+1
    tl=tvec(ll-1)+hvec(ll-1)*c;
    %integral along circular contour
    int=0;
    for mm=1:length(Z)
        U(:,:,mm)=transpose((I-hvec(ll-2)*A*Z(mm))\transpose(((repmat(U(:,end,mm),1,s)+(C(:,:,ll-1)-repmat(C(:,end,ll-2),1,s))))));
        int=int+weights(mm)*(Z(mm)^(alp-1))*transpose((I-hvec(ll-1)*Z(mm)*A)\transpose(repmat(U(:,end,mm),1,s)));
    end
    P1=((rf(f,tl,J,s)-int)*transpose((hvec(ll-1)*A)^alp)+repmat(C(:,end,ll-1),1,s))*transpose(invV);
    P2=idstf(P1,J,s);
    P3=1./(1-transpose(eig(A).^alp).*(hvec(ll-1)^alp*reshape(D,(J-1)^2,1))).*P2;
    C(:,:,ll)=dstf(P3,J,s)*transpose(V);
    E =[E norm(reshape(C(:,end,ll),J-1,J-1)-sol(tvec(ll)),'fro')*h];
end

if n0<N
    
    %          dt=c(1)*dt;
    [X,W,Q1,Q2s]=quadrature_cqw_grad(1-alp,tolquad,n0,tvec(end-n0),dt,hmax,B,RK);
    X=X';
    W=W';
    %Correction factor for memory term
    factor=ones(1,length(X));
    RVn=zeros(s,s*length(X));
    rn=ones(1,length(X));
    for ll=1:n0
        %         for mm=1:length(X)
        %             factor(mm)=factor(mm)*R(-hvec(ll)*X(mm));
        %         end
        factor=factor./R(-hvec(ll)*X);
        
    end
    %For ODEs compression of the memory
    Y=zeros((J-1)^2,s,length(X));
end

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
    U=zeros((J-1)^2,s,length(Z));
    for ll=1:n0-1
        %integral along circular contour
        int=0;
        for mm=1:length(Z)
            U(:,:,mm)=transpose((I-hloc(ll)*A*Z(mm))\transpose(repmat(U(:,end,mm),1,s)+C(:,:,n-n0+ll+1)-repmat(C(:,end,n-n0+ll),1,s)));
            if ll==n0-1
                int = int + weights(mm)*(Z(mm)^(alp-1))*transpose((I-hvec(n)*Z(mm)*A)\transpose(repmat(U(:,end,mm),1,s)));
            end
        end
        %         local=int;
    end
    
    
    local=int;
    memory=0;
    %Memory Term for n>n0
    %     for mm=1:length(X)
    %         factor(mm)=factor(mm)*(R(-hvec(n)*X(mm))/R(-hvec(n-n0)*X(mm)));%？？？？
    %     end
    factor=factor.*(R(-hvec(n-n0)*X)./R(-hvec(n)*X));
    %Update ODEs up to t_{n-n0} and quadrature
    
    for mm=1:length(X)
        
        Y(:,:,mm)=transpose((I+hvec(n-n0)*X(mm)*A)\transpose(repmat(Y(:,end,mm),1,s)+C(:,:,n-n0+1)-repmat(C(:,end,n-n0),1,s)));
        RVn(:,(mm-1)*s+1:mm*s)=(I+hvec(n)*X(mm)*A)\Mat1;
        rn(mm)=R(-hvec(n)*X(mm));
        memory=memory+transpose((RVn(:,(mm-1)*s+1:mm*s)/rn(mm))*transpose((W(mm)/factor(mm)*Y(:,:,mm))));
    end
    %     memory=memory+W.*factor*Y(:,mm);
    tn=tvec(n)+hvec(n)*c;
    P1=((rf(f,tn,J,s)-(local+memory))*transpose((hvec(n)*A)^alp)+repmat(C(:,end,n),1,s))*transpose(invV);
    P2=idstf(P1,J,s);
    P3=1./(1-transpose(eig(A)).^alp.*(hvec(n)^alp*reshape(D,(J-1)^2,1))).*P2;
    C(:,:,n+1)=real(dstf(P3,J,s)*transpose(V));
    E =[E norm(reshape(C(:,end,n+1),J-1,J-1)-sol(tvec(n+1)),'fro')*h];
end
%%

%Error

[eamax]=max(E);
eTf=E(end);




end
