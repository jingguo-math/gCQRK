% Test the convergence rate of gCQ RK for K_a(\partial_t)f, where K(z)=1/(z^alp+r)
clc
clf
close all
clear

% Add path to quadrature functions (relative path from current location)
addpath('../../quadratures');  

% Numerical method parameters
RK = 2;         % Runge-Kutta method order
Ord = 3;         % Expected convergence order
Tf = 1;          % Final time (original value was 1)

% Problem parameters
alp_values = 0.3:0.1:0.7; % Vector of alpha values to test
bet =0.2;       % Beta parameter for the test function
Nvec = 2*2.^[2:9]; % Vector of discretization points (N values)

% Test function f(t) = t^bet
f = @(t) t.^bet;

% Initialize cell array to store errors for each alpha value
all_errors = cell(length(alp_values), 1);

% Main loop over different alpha values
for alp_idx = 1:length(alp_values)
    alp = alp_values(alp_idx);
    
    % Kernel function 
    Kfun = @(z)((z^alp+speye(size(z)))\speye(size(z)));
    
    % Exact solution (using Mittag-Leffler function)
    sol = @(t)gamma(bet+1)*t.^(alp+bet).*ml(-t.^alp,alp,alp+bet+1);
    
    % Gradient parameter (adapts to problem parameters)
    grad = max(1, Ord/(alp+bet));
    
    % Initialize error vector
    E = [];
    
    % Convergence study loop over different N values
    MaxIt = length(Nvec);
    for j = 1:MaxIt
        N = Nvec(j);
        % Call the generalized composite quadrature Runge-Kutta method
        [U,e] = cqrk_varn0_genkTrap(sol,f,Kfun,RK,N,Tf,alp,grad);
        E = [E,e]; % Store errors
    end
    
    % Store errors for this alpha value
    all_errors{alp_idx} = E;
    
    % Compute convergence rates (log2 of error ratios)
    rate = log2(E(1:end-1)./E(2:end));
    
    %% Display results for current alpha
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

% Define line styles and colors for different alpha values
line_styles = {'d-', 'p:', '*--', 's-.', 'o--'};
colors = [
    0 0.4470 0.7410;    % Blue
    0.8500 0.3250 0.0980; % Orange
    0.9290 0.6940 0.1250; % Yellow
    0.4940 0.1840 0.5560; % Purple
    0.4660 0.6740 0.1880  % Green
    ];

% Plot reference line for expected convergence rate (slope = -Ord)
loglog(Nvec, all_errors{2}(3)*(Nvec(3)./Nvec).^Ord, 'k', 'LineWidth', 2);
hold on

% Create legend entries (reference line + alpha values)
legend_entries = cell(1, length(alp_values) + 1);
legend_entries{1} = sprintf('$slope=%.0f$', -Ord); % Reference slope

% Plot results for each alpha value
for alp_idx = 1:length(alp_values)
    legend_entries{alp_idx+1} = sprintf('$\\alpha=%.1f$', alp_values(alp_idx));
    loglog(Nvec, all_errors{alp_idx}, line_styles{alp_idx}, ...
        'LineWidth', 2, 'MarkerSize', 12, 'Color', colors(alp_idx,:));
    hold on
end

% Set axis labels (LaTeX formatted)
xlabel('$N$', 'FontSize', 30, 'Interpreter', 'Latex');
ylabel('Maximum Absolute Error', 'FontSize', 30, 'Interpreter', 'Latex');

% Set legend (positioned in southwest corner)
legend(legend_entries, 'Location', 'southwest', 'FontSize', 23, 'Interpreter', 'Latex');

% Set x-axis ticks and labels
xticks(Nvec);
set(gca, 'XTickLabel', Nvec, 'FontName', 'Times', 'FontSize', 26);

% Set y-axis ticks in scientific notation
yticks(10.^(-10:2:-2));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('10^{%d}', x), -10:2:-2, ...
    'UniformOutput', false), 'FontName', 'Times', 'FontSize', 26);

% Set axis limits
set(gca, 'XLim', [Nvec(1)*0.9, Nvec(end)*1.1]);

% Set y-axis limits based on actual data range
all_errors_vector = cell2mat(all_errors);
valid_errors = all_errors_vector(all_errors_vector > 0);
if ~isempty(valid_errors)
    y_lower = max(min(valid_errors)*0.5, 1e-15);
    y_upper = min(max(valid_errors)*2, 1e10);
    set(gca, 'YLim', [y_lower, y_upper]);
end

% Set figure size and position
set(gcf, 'Position', [100, 100, 700, 500]);

% Customize grid appearance
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';
ax.GridColor = [0.7, 0.7, 0.7];
ax.GridAlpha = 0.7;
ax.LineWidth = 1.5;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

% Final plot formatting
box on;
hold off;