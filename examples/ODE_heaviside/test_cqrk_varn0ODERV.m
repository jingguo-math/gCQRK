% approximate D_t^\alp u+mu*u=f(f is a given function) at prescribed general time mesh tn, n=1,...,N
% the exact solution has two different singularities
% Clear command window, variables, figures
clc
clear
close all
clf

% Add path to quadrature functions (relative path from current location)
addpath('../../quadratures');  
% Numerical method parameters
RK=2;       % Runge-Kutta method order (2nd order)
Ord=3;      % Expected convergence order
Tf=1;       % Final simulation time

% Problem parameters
alp_values = [0.25 0.5 0.75]; % Fractional derivative orders to test
bet1=0.5;   % First singularity exponent
bet2=0.9;   % Second singularity exponent
grad1=Ord/bet1; % Grading parameter for first singularity
grad2=Ord/bet2; % Grading parameter for second singularity

mu=1;       % Coefficient in the differential equation
Nvec=[8 16 32 64 128 256 512 1024 2048]; % Grid sizes to test
sigm=0.72;  % Time where second singularity occurs

% Heaviside function definition
H=@(t)t>=0;

% Exact solution with two singularities:
% - First term (1 + t^β1) for t ∈ [0,Tf]
% - Second term (t-σm)^β2 activates at t = σm
bet=bet1; % Current beta value being tested
u = @(t)1+t.^bet1+H(t-sigm).*(t-sigm).^bet2;

% Initialize storage for errors (one cell per alpha value)
all_errors = cell(length(alp_values), 1);

% Main loop over different fractional orders (α values)
for alp_idx = 1:length(alp_values)
    alp = alp_values(alp_idx);
    
    % Right-hand side function f(t) for the fractional ODE:
    % f(t) = μu(t) + D_t^α u(t) (fractional derivative of exact solution)
    f=@(t) mu*u(t)+gamma(bet1+1)*t.^(bet1-alp)/gamma(bet1-alp+1)+...
        gamma(bet2+1)*H(t-sigm).*(t-sigm).^(bet2-alp)/gamma(bet2-alp+1);
    
    E=[]; % Initialize error vector for current alpha
    
    % Convergence study loop over different grid sizes
    for N=Nvec
        % Call the fractional ODE solver:
        % - u: exact solution
        % - f: right-hand side
        % - RK: Runge-Kutta order
        % - N: number of time steps
        % - Tf: final time
        % - alp: fractional order
        % - mu: equation coefficient
        % - sigm: singularity location
        % - grad1, grad2: grading parameters
        [C,eamax]=cqrk_varn0ODERV(u,f,RK,N,Tf,alp,mu,sigm,grad1,grad2);
        E=[E;eamax]; % Store maximum absolute error
    end
    
    % Store errors for this alpha value
    all_errors{alp_idx} = E;
    
    % Compute convergence rates (log2 of error ratios)
    rate = log2(E(1:end-1)./E(2:end));
    
    %% Display convergence results for current alpha
    fprintf('\nResults for alp = %.1f:\n', bet);
    disp('RK method:');
    disp(RK);
    disp('alp:');
    disp(alp);
    disp('bet:');
    disp(bet);
    disp('Convergence rates:');
    disp(rate);
end

%% Plot convergence results
figure(1);

% Line styles and colors for different alpha values
line_styles = {'d-', 'p:', '*--', 's-.', 'o--'};
colors = [
    0 0.4470 0.7410;    % Blue
    0.8500 0.3250 0.0980; % Orange
    0.9290 0.6940 0.1250; % Yellow
    0.4940 0.1840 0.5560; % Purple
    0.4660 0.6740 0.1880  % Green
    ];

% Plot reference line showing expected convergence rate (slope = -Ord)
% Uses 6th error point as reference for scaling
loglog(Nvec, all_errors{end}(6)*(Nvec(6)./Nvec).^Ord, 'k', 'LineWidth', 2);
hold on

% Create legend entries (reference line + alpha values)
legend_entries = cell(1, length(alp_values) + 1);
legend_entries{1} = sprintf('$slope=%.0f$', -Ord); % Reference slope

% Plot error curves for each alpha value
for alp_idx = 1:length(alp_values)
    legend_entries{alp_idx+1} = sprintf('$\\alpha=%.2f$', alp_values(alp_idx));
    loglog(Nvec, all_errors{alp_idx}, line_styles{alp_idx}, ...
        'LineWidth', 2, 'MarkerSize', 12, 'Color', colors(alp_idx,:));
    hold on
end

% Figure formatting
xlabel('$N$', 'FontSize', 30, 'Interpreter', 'Latex');
ylabel('Maximum Absolute Error', 'FontSize', 30, 'Interpreter', 'Latex');
legend(legend_entries, 'Location', 'southwest', 'FontSize', 23, 'Interpreter', 'Latex');

% Axis ticks and labels
xticks(Nvec);
set(gca, 'XTickLabel', Nvec, 'FontName', 'Times', 'FontSize', 26);
yticks(10.^(-10:2:-2));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('10^{%d}', x), -10:2:-2, ...
    'UniformOutput', false), 'FontName', 'Times', 'FontSize', 26);

% Axis limits
set(gca, 'XLim', [Nvec(1)*0.9, Nvec(end)*1.1]);
all_errors_vector = cell2mat(all_errors);
valid_errors = all_errors_vector(all_errors_vector > 0);
if ~isempty(valid_errors)
    y_lower = max(min(valid_errors)*0.5, 1e-15);
    y_upper = min(max(valid_errors)*10, 1e10);
    set(gca, 'YLim', [y_lower, y_upper]);
end

% Figure size and grid formatting
set(gcf, 'Position', [100, 100, 700, 500]);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';
ax.GridColor = [0.7, 0.7, 0.7];
ax.GridAlpha = 0.7;
ax.LineWidth = 1.5;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

box on;
hold off;