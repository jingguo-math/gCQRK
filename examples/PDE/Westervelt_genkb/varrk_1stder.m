function [f] = varrk_1stder(dt,A,s,U,J)
%Compute the variable step Runge-Kutta for the first order derivative \partial_t U at
%t_n
%Input: dt:time step size; A:coeffient matrix in the Butcher's table;
%U:vector [u^{n-1} u^n],u^{n-1} u^n is a column verctor with (s*J)
%components
%Output: gCQ Runge-Kutta approximation of the first order derivative \partial_t U 
%==============Written by Jing Guo, 28.05.2024=============================
rhsV=U(:,end)-repmat(U((s-1)*J+1:s*J,end-1),s,1);
f=transpose((dt*A)\transpose(reshape(rhsV,J,s)));
end

