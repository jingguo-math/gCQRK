function [C,eamax]=cqrk_varn0ODERV(sol,f,RK,N,Tf,alpha,mu,r,grad1,grad2)

% approximate D_t^\alpha u+mu*u=f(f is a given function) at prescribed general time mesh tn, n=1,...,N
%using quadrature on the hyperbola to compute the history part
%It uses gCQ based on a Runge--Kutta method
%ref1: Generalized convolution quadrature with variable time
% stepping. Part II: Algorithm and numerical results
%=============== Modified by Jing Guo, May, 2024===============
%%
tolquad=1e-14;%tolerance of the quadrature rule
[s,A,c,b,eiginvA,V,invV]=setsolver(RK); %Buther tableau and eigenvalue decomposition of A^(-1)
Mat1=zeros(s);
Mat1(:,end)=ones(s,1);
I=eye(s); 
if RK == 1
    R = @(z) 1./(1 - z);
elseif RK == 2
    R = @(z) (2*z + 6)./(z.^2 - 4*z + 6);
else
    R = @(z) 1 + z*b'*((eye(s) - z*A)\ones(s,1)); % stability function
end

if grad1>40
    n0=N;
else
    n0=min(5,N); %Direct steps for RK-based gCQ
end
B = 2;
%% construct the time mesh
tvec=zeros(1,N+1);
N1=floor(N*r/Tf);
step_l1= (1:1:N1).^(grad1-1).*(N1:-1:1).^(grad2-1)*N1^(-grad1-grad2+1);
hvec_l1=r/sum(step_l1)*step_l1;
for l=1:N1
    tvec(l+1)=tvec(l)+hvec_l1(l);
end
tvec(N1+2:N+1)=tvec(N1+1)+((1:N-N1)./(N-N1)).^(grad2)*(Tf-r);
hvec = tvec(2:end)-tvec(1:end-1); %temporal step
%% Figure of the mesh
figure(1); clf%time stepsize figure
semilogy(tvec(2:end),hvec,'.-','Linewidth',1.5,'Markersize',15);
xlim([0 1.001]);
ylim([min(hvec) max(hvec)*1.2]);
xticks([0 0.2 0.4 0.6 0.8 1]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',26)
hp = xlabel('$t$');
set(hp,'FontSize',30,'Interpreter','Latex')
hp = ylabel('$\tau$');
set(hp,'FontSize',30,'Interpreter','Latex')


dt=min(hvec(1:end-n0));%minimum time step in the history part
hmax=max(hvec(1:end-n0));%maximum time step in the history part
%Approximation of convolution at first n0 steps
C=zeros(s,N+1);
% C(:,1)=sol(0)*ones(s,1);
Cin=sol(0)*ones(s,1);
%Parameters for elliptic quadrature
if RK==1
    circ_radi=3;%choice of the radius of the circle contour
else
    circ_radi=1;
end
qvec=zeros(1,N-n0);
tn=hvec(1)*c;
C(:,2) = (I+mu*(hvec(1)*A)^alpha)\((hvec(1)*A)^alpha*(f(tn)-mu*Cin)+(C(end,1)*ones(s,1)));
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
U=zeros(s,length(Z));
for ll=3:n0+1
    tl=tvec(ll-1)+hvec(ll-1)*c;
    %integral along circular contour
    int=0;
    for mm=1:length(Z)
        U(:,mm)=(I-hvec(ll-2)*A*Z(mm))\(U(end,mm)*ones(s,1)+(C(:,ll-1)-C(end,ll-2)*ones(s,1)));
        int=int+weights(mm)*(Z(mm)^(alpha-1))*((I-hvec(ll-1)*Z(mm)*A)\(U(end,mm)*ones(s,1)));
    end
    C(:,ll)=(I+mu*(hvec(ll-1)*A)^alpha)\((hvec(ll-1)*A)^alpha*f(tl)+...
        (C(end,ll-1)*ones(s,1))-mu*(hvec(ll-1)*A)^alpha*Cin-(hvec(ll-1)*A)^alpha*int);        
end
%%
%Approximation of convolution for n>n0

%Quadrature for the memory term
if n0<N
    
    %          dt=c(1)*dt;
    [X,W,Q1,Q2s]=quadrature_cqw_grad(1-alpha,tolquad,n0,tvec(end-n0),dt,hmax,B,RK);
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
    Y=zeros(s,length(X));
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
 
    U=zeros(s,length(Z));
    for ll=1:n0-1
        %integral along circular contour
        for mm=1:length(Z)
            U(:,mm)=(I-hloc(ll)*A*Z(mm))^(-1)*(U(end,mm)*ones(s,1)+C(:,n-n0+ll+1)-C(end,n-n0+ll)*ones(s,1));
        end
    end
    int=zeros(s,1);
    for mm=1:length(Z)
        int = int + weights(mm)*(Z(mm)^(alpha-1))*((I-hvec(n)*Z(mm)*A)\(U(end,mm)*ones(s,1)));
    end
    local=int;
    memory=zeros(s,1);
    factor=factor.*(R(-hvec(n-n0)*X)./R(-hvec(n)*X));
    %Update ODEs up to t_{n-n0} and quadrature
    
    for mm=1:length(X)
        
        Y(:,mm)=(I+hvec(n-n0)*X(mm)*A)\(ones(s,1)*Y(end,mm)+C(:,n-n0+1)-C(end,n-n0)*ones(s,1));
        RVn(:,(mm-1)*s+1:mm*s)=(I+hvec(n)*X(mm)*A)\Mat1;
        rn(mm)=R(-hvec(n)*X(mm));
        memory=memory+RVn(:,(mm-1)*s+1:mm*s)/rn(mm)*(W(mm)/factor(mm)*Y(:,mm));
    end
    tn=tvec(n)+hvec(n)*c;
    C(:,n+1)= ((hvec(n)*A)^-alpha+mu*I)\(f(tn)-mu*Cin-(local+memory)+(hvec(n)*A)^-alpha*(C(end,n)*ones(s,1)));

end
C=real(C+Cin);

%Error

solvec = sol(tvec(1:end));
err = abs(C(end,:)-solvec);

[eamax,indea]=max(err);
figure(3); clf;
plot(tvec,err,'*-');


figure(20); clf;
plot(tvec(1:end),C(end,:),'b.',tvec(1:end),solvec,'r-');




end
