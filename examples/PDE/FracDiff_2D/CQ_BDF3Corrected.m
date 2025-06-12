function [C,eamax,etf]=CQ_BDF3Corrected(sol,f,N,Tf,J,alp,h,c_u)

% approximate the 2d PDE D_t^\alpha u-c_u*Lap*u=f(f is a given function) at prescribed general time mesh tn, n=1,...,N
%using quadrature on the hyperbola to compute the history part

%It uses corrected CQ  based on the BDF3
%ref1:CORRECTION OF HIGH-ORDER BDF CONVOLUTION QUADRATURE FOR FRACTIONAL
%EVOLUTION EQUATIONS (7)
%ref2: Fast High-Order Compact Exponential Time Differencing Runge–Kutta Methods for Second-Order Semilinear Parabolic Equations
%=============== Modified by Jing Guo, May, 2024===============

%%

E=[];
Ks = @(z) (z).^alp;
dlt= @(z)(1-z)+1/2*(1-z).^2+1/3*(1-z).^3;
dt=Tf/N;
Lambda=eps^(1/2/N);
Lam=Lambda.^(0:N);
L=N+1;
w=fft(Ks(dlt(Lambda*exp(1i*2*pi*(0:N)./L))/dt))./(L*Lam);
tvec = ((0:N)*Tf/N);
% dt = tvec(2)-tvec(1);
D=-4/h^2*(sin((1:(J-1))'*pi/(2*J))).^2;
D=c_u*(D+D'+1/12*(h^2+h^2)*D*D')./(1+1/12*(h^2*D+h^2*D'));%ref 2: Page 6
C=zeros(J-1,J-1,N+1);
C(:,:,1)=sol(0);


MatE=1./(w(1)-D);
idstx=@(x)transpose(idst(transpose(idst(x))));
dstx=@(x)transpose(dst(transpose(dst(x))));
Dstx=@(x)real(dstx(MatE.*idstx(x)));
Lx=@(x)(8*x+[zeros(J-1,1)  x(:,1:end-1)]+[x(:,2:end) zeros(J-1,1)]+...
        [zeros(1,J-1);  x(1:end-1,:)]+[x(2:end,:); zeros(1, J-1)])/12;
     rhs=(11/12*(Lx(C(:,:,1))+f(0))+f(tvec(2))+1/24*(f(tvec(1)+dt)-f(tvec(1)-dt)))+w(1)*C(:,:,1);
%   rhs=(11/12*(Lx(C(:,:,1))+f(0))+f(tvec(2))+1/12*(f(tvec(2))))+w(1)*C(:,:,1);
% rhs=(11/12*(Lx(C(:,:,1))+f(0))+f(tvec(2))+1/12*(f(tvec(2))-(tvec(1))))+w(1)*C(:,:,1);
%  rhs=(11/12*(Lx(C(:,:,1))+f(0))+f(tvec(2)))+w(1)*C(:,:,1);
C(:,:,2)=Dstx(rhs);
E =[E norm(C(:,:,2)-sol(tvec(2)),'fro')*h];
rhs=(-5/12*(Lx(C(:,:,1))+f(0))+f(tvec(3))-w(2)*(C(:,:,2)-C(:,:,1)))+w(1)*C(:,:,1);
C(:,:,3)=Dstx(rhs);
E =[E norm(C(:,:,3)-sol(tvec(3)),'fro')*h];



for ll=4:N+1
    sum1=0;
    for k=2:ll-1
    sum1=sum1+(C(:,:,k)-C(:,:,1))*w(ll-k+1);
    end
   rhs=f(tvec(ll))-sum1+w(1)*C(:,:,1);
    C(:,:,ll)=Dstx(rhs);
   E =[E norm(C(:,:,ll)-sol(tvec(ll)),'fro')*h];
end


[eamax]=max(E);
etf=E(end);



end
