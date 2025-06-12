function [Z,weights]=elliptic_quadrature(minpole,maxpole,nq)
% Written by Maria Lopez-Fernandez 2022.

m=minpole; 
M=maxpole;
shif=m/10;

if  M/m < 1.1
    %Simple trapezoidal quadrature if M/m is close to 1
    mid=M; r=M;
    gam=@(t)mid+r*exp(1i*t);
 gam=@(t) gam(t)+shif;
    dgam=@(t) r*exp(1i*t);
    tnodes=-pi:2*pi/nq:pi;
    tnodes=tnodes(1:end-1);
    Z = gam(tnodes);   
    weights = -dgam(tnodes)/(nq);
%     figure(1)
%     hold on
%     plot(Z,'.-')
%     axis equal
%      figure(2)    
%     hold on
%     plot(weights,'.-')    
%     axis equal

else
    kappa = M*m/(M-m);
    rho = (M+kappa)/(m+kappa);
    cq = (m+kappa)*sqrt(rho);
    k = (sqrt(rho)-1)/(sqrt(rho)+1);
    L = -log(k)/pi;
    [K,Kp] = ellipkkp(L);

    %Nodes and weights
    t = .5i*Kp - K + (.5:nq)*4*K/nq;
    [u cn dn] = ellipjc(t,L);
    Z = cq*((1/k+u)./(1/k-u))-kappa;
%     Z = Z+1;
    Z = Z+shif;
    weights = cn.*dn./(1/k-u).^2;
    weights = 4*K*cq*weights/(pi*1i*k*nq);
%     figure(1)    
%     hold on
%     plot(Z,'.-')
%     axis equal
%     figure(2)    
%     hold on
%     plot(weights,'.-')    
%     axis equal
end

