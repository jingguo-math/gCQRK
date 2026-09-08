function [U, E, NQ] = cqrk_genkTrap(sol, f, Kfun, RK, N, Tf, alpha, grad, bet)
%Fast RK-gCQ for Example 2 with K_a(z)=1/(z^alpha+1).
%   The local-memory splitting uses n0=1. The history integral is
%   discretized after x=exp(mu/alpha) by the trapezoidal rule prescribed
%   in Example 2.
%
%   [U,E,NQ] returns the RK stage values, the maximum nodal error, and the
%   number NQ=2*Mtilde+1 of  quadrature nodes. 

%=============== Modified by Jing Guo, 08-09-2026===============

validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(Tf, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(alpha, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<', 1});
validateattributes(grad, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(bet, {'numeric'}, {'scalar', 'real', 'finite'});
if alpha + bet <= 0
    error('cqrk_varn0_genkTrapV3:InvalidRegularity', ...
        'Example 2 requires alpha+bet>0.');
end

[s, A, c, ~, eiginvA, V, invV] = setsolver(RK);
I = eye(s);
one = ones(s, 1);

% Graded mesh t_n=T(n/N)^grad.
tvec = Tf * ((0:N) / N).^grad;
hvec = diff(tvec);
U = zeros(s, N + 1);
err = zeros(1, N);
NQ = 0;

% For n0=1 the entire local contribution is the current-step term
% K((h_n A)^(-1)) f_n; no local contour quadrature is required.
K0 = V * diag(Kfun(eiginvA / hvec(1))) * invV;
fstage = f(tvec(1) + c * hvec(1));
U(:, 2) = real(K0 * fstage);
err(1) = U(end, 2) - sol(tvec(2));

if N >= 2
    % Example 2: d=0.9*min(pi(1-alpha),pi*alpha/2) and
    % Mtilde=ceil(((r*log(N)+2)^2)/(2*pi*d)), where
    % r=min(p,q+1+alpha,grad*(alpha+bet)).
    if (isnumeric(RK) && any(RK == [1, 2, 3])) || ...
            (ischar(RK) && strncmp(RK, 'RadauIIA', 8)) || ...
            (isstring(RK) && startsWith(RK, "RadauIIA"))
        p = 2 * s - 1;
        q = s;
    else
        p = 2 * (s - 1);
        q = s - 1;
    end
    d = 0.9 * min(pi * (1 - alpha), pi * alpha / 2);
    r = min([p, q + 1 + alpha, grad * (alpha + bet)]);
    Mtilde = max(1, ceil((r * log(N) + 2)^2 / (2 * pi * d)));
    htrap = sqrt(2 * pi * d / Mtilde);
    [X, W] = Trap_quadrature(alpha, Mtilde, Mtilde, htrap);
    X = X(:).';
    W = W(:).';
    NQ = numel(X);

    % Q_j^his(x), initialized by Q_0^his=0. At step n, update it only
    % through j=n-1 and then apply Rvec(-h_n*x)e_s^T, as n0=1.
    Qhis = zeros(s, NQ);
    for n = 2:N
        hprev = hvec(n - 1);
        fprev = f(tvec(n - 1) + c * hprev);
        memory = zeros(s, 1);

        for ell = 1:NQ
            Qhis(:, ell) = (I + hprev * X(ell) * A) \ ...
                (one * Qhis(end, ell) + hprev * A * fprev);

            Rvec = (I + hvec(n) * X(ell) * A) \ one;
            memory = memory + W(ell) * Rvec * Qhis(end, ell);
        end

        fstage = f(tvec(n) + c * hvec(n));
        K0 = V * diag(Kfun(eiginvA / hvec(n))) * invV;
        U(:, n + 1) = real(K0 * fstage + memory);
        err(n) = U(end, n + 1) - sol(tvec(n + 1));
    end
end

if any(~isfinite(U(:))) || any(~isfinite(err))
    error('cqrk_varn0_genkTrapV3:NonfiniteResult', ...
        'Computed solution or error contains NaN or Inf.');
end
E = max(abs(err));
end
