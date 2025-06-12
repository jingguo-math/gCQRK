function [X,weights]=Trap_quadrature(alpha,M,N,tau)
% trapexoidal quadrature for
% z^{\alpha-1}(z^\alpha+A)^{-1}
n=-M:N;
X=exp(n*tau/alpha);
weights=sin(pi*alpha)/pi*tau/alpha.*X./(exp(n*tau)+2*cos(pi*alpha)+exp(-n*tau));
end

