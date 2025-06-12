function [X,W,Q1,Q2s]=quadrature_cqw_grad(alpha,tol,n0,T,dt,hmax,B,RK)
Qmax = 200;
% weights 0, \dots, n_0 computed directly and n_0+1 etc use quadrature
%adjust the error bound by replacing some dt by hmax

%constants for estimates
if (RK == 1) %BDF1 
    cc = 1;x0 = 1;
    b = .5; gam = -log(1-b)/b; Cq = exp(gam*b); Cq1 = 1;
    rn = @(z) 1./abs(1-z); qn = @(z) 1./abs(1-z);
else %currently just 2-stage Radau
    cc = 1/3.4641 ;x0 = 2*cc;
    b = 1; gam = 1.0735; Cq = 1.4733; 
    Cq1 = 1;%need to check what Cq1 is for Radau
    rn = @(z) abs((2*z+6)./(z.^2-4*z+6)); qn = @(z) sqrt((9/2)^2+abs(3/2-z).^2)./abs(z.^2-4*z+6); 
end

A = 0;
err = inf;
while (err > tol/3)%p15 lem. 14
    A = A+.125;
% err = (dt^alpha/(gamma(alpha)*gamma(1-alpha)))*integral(@(x) qn(-x).*(rn(-x).^(n0+1)).*x.^(-alpha),A,inf);   
 err = (dt^(alpha-1)*hmax/(n0*gamma(alpha)*gamma(1-alpha)))*A^(-alpha)*(x0+cc*A)^(-n0);
end
%Atmp = (dt^alpha/(.333*tol*n0*gamma(1-alpha)*cc^(n0+2)))^(1/(n0+1+alpha)); A = max(1,A);
%disp([A Atmp])    

L = A/dt; % to achieve a truncation error bounded by tol

% Do Gauss-Jacobi on interval [0,L0]
L0 = 4/T;
J = floor(log(L/L0)/log(1+B)); B = (L/L0)^(1/J)-1; %improves things a little bit
%number of subintervals
%L = L0*(1+B)^J; 
rho_max = 1+2*b/(L0*hmax)+sqrt((2*b/(L0*hmax))^2+4*b/(L0*hmax));
%  rho_max = 1+2*b/(L0*dt)+sqrt((2*b/(L0*dt))^2+4*b/(L0*dt));
err = inf;
Q1 = 0;

while (Q1 < Qmax && err > tol/3)%p16 thm 17
    Q1 = Q1+1;
    rho_opt = 4*Q1/(gam*T*L0) + sqrt(1+(4*Q1/(gam*T*L0))^2);
    if (rho_opt < rho_max)
%         disp('rho_opt < rho_max')
%         disp([rho_opt rho_max])
        err = Cq*(hmax/((1-alpha)*gamma(alpha)*gamma(1-alpha)))*(L0^(1-alpha))*(1+gam*T*L0/(4*Q1))*(exp(1)*gam*T*L0/(8*Q1))^(2*Q1);
    else
        err = Cq*L0^(1-alpha)*(hmax/((1-alpha)*gamma(alpha)*gamma(1-alpha)))*((exp(gam*T*b/hmax)*rho_max^(-2*Q1+1)/(rho_max-1)));
%         disp('rho_opt >= rho_max')
%         disp([rho_opt rho_max])
    end
end


g_epss = @(epss) 1+(2/B)*(1-epss)+sqrt((1+(2/B)*(1-epss)).^2-1);
er_epss = @(epss,Lj,Q) (min(Cq1,(x0+cc*Lj*dt*epss).^(-n0-2))).*(epss.^(-alpha)).*(g_epss(epss).^(-2*Q+1))./(g_epss(epss)-1);
%er_epss = @(epss,Lj,Q) Cq*(epss.^(-alpha)).*(g_epss(epss).^(-2*Q+1))./(g_epss(epss)-1);
Q2s = zeros(J,1);
%Computation of all quadrature weights and nodes 

%Nodes and weights for the integral in [0,L0]
[xjac,wjac] = gaussj(Q1,0,-alpha); 
xjac = (xjac+1)*L0/2; 
wjac = wjac*(2/L0).^(-1+alpha);
X=xjac; W=wjac;

%Nodes and weights in [L_{j-1},L_j], j=1,...,J

for j = 0:(J-1)% p18 thm 19
    Lj = L0*((1+B)^j);
   
    Q2 = 0;
    err = inf;
    while (Q2 < Qmax && err > tol/3)
        Q2 = Q2+1;
        err = min(er_epss((1:99)/100,Lj,Q2));                
        err = hmax*(1/(2*gamma(alpha)*gamma(1-alpha)))*B*(Lj^(1-alpha))*err;    
    end
    Q2s(j+1) = Q2;
    [xg,wg] = gauss(Q2); 
    
    xgj = (xg+1)*B*Lj/2+Lj;
    wgj = wg.'*B*Lj/2;
    wgj=wgj.*xgj.^(-alpha);
    X=[X;xgj];
    W=[W;wgj];
end
W=W/(gamma(1-alpha)*gamma(alpha ));