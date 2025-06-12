% test the convergnce raten of gCQ RK for K_a(\partial_t)f, where K(z)=1/(z+1)^alp 
clc
clf
close all
clear

% Add path to quadrature functions (relative path from current location)
addpath('../../quadratures'); 

% Numerical method parameters
RK = 2;         % Runge-Kutta method order
Ord = 3;        % Expected convergence order
Tf = 1;         % Final simulation time

% Problem parameters
alp_values = [0.3 0.4 0.5 0.6 0.7]; % Vector of α values to test
bet = 0.2;      % Exponent for the test function f(t) = t^β
Nvec = 2*2.^[2:9]; % Vector of discretization points (N values)

% Test function definition f(t) = t^β
f = @(t) t.^bet;

% Initialize cell array to store errors for each α value
all_errors = cell(length(alp_values), 1);

% Main convergence study loop over different α values
for alp_idx = 1:length(alp_values)
    alp = alp_values(alp_idx);
    
    % Kernel function K(z) = (z+I)^(-α) (discrete version)
    Kfun = @(z)((z+speye(size(z)))^alp\speye(size(z)));
    
    % Gradient parameter adapts to problem parameters
    grad = max(1, Ord/(alp+bet));
    
    % Initialize error vector
    E = [];
    
    % Reference solution with finest resolution
    U = cqrk_varn0_genkb(f, Kfun, RK, Nvec(1), Tf, alp, grad);
    
    % Convergence study loop over different N values
    MaxIt = length(Nvec);
    for j = 1:MaxIt
        N = Nvec(j)*2;  % Double resolution for reference solution
        
        % Compute reference solution
        U_ref = cqrk_varn0_genkb(f, Kfun, RK, N, Tf, alp, grad);
        
        % Compute maximum absolute error between current and reference solutions
        % (compares every other point for matching grid sizes)
        E = [E, max(abs(U_ref(end,1:2:end) - U(end,:))];
        
        % Update reference solution for next iteration
        U = U_ref;
    end
    
    % Store errors for this α value
    all_errors{alp_idx} = E;
    
    % Compute convergence rates (log2 of error ratios between successive N)
    rate = log2(E(1:end-1)./E(2:end));
    
    %% Display results for current α value
    disp('RK method:');
    disp(RK);
    disp('α:');
    disp(alp);
    disp('β:');
    disp(bet);
    disp('Convergence rates:');
    disp(rate);
end

%% Plot convergence results
figure(1);

% Define line styles and marker types for different α values
line_styles = {'d-', 'p:', '*--', 's-.', 'o--'}; % diamond, pentagon, asterisk, square, circle

% Define color palette for plots (MATLAB default colors)
colors = [
    0 0.4470 0.7410;    % Blue
    0.8500 0.3250 0.0980; % Orange
    0.9290 0.6940 0.1250; % Yellow
    0.4940 0.1840 0.5560; % Purple
    0.4660 0.6740 0.1880  % Green
    ];

% Plot reference line showing expected convergence rate (slope = -Ord)
% Uses 5th error point as reference for scaling
loglog(Nvec, all_errors{end}(5)*(Nvec(5)./Nvec).^Ord, 'k', 'LineWidth', 2);
hold on

% Create legend entries (reference line + α values)
legend_entries = cell(1, length(alp_values) + 1);
legend_entries{1} = sprintf('$slope=%.0f$', -Ord); % Reference slope entry

% Plot error curves for each α value
for alp_idx = 1:length(alp_values)
    legend_entries{alp_idx+1} = sprintf('$\\alpha=%.1f$', alp_values(alp_idx));
    loglog(Nvec, all_errors{alp_idx}, line_styles{alp_idx}, ...
        'LineWidth', 2, 'MarkerSize', 12, 'Color', colors(alp_idx,:));
    hold on
end

% Set axis labels with LaTeX formatting
xlabel('$N$ (Number of time steps)', 'FontSize', 30, 'Interpreter', 'Latex');
ylabel('Maximum Absolute Error', 'FontSize', 30, 'Interpreter', 'Latex');

% Configure legend (positioned in southwest corner)
legend(legend_entries, 'Location', 'southwest', 'FontSize', 23, 'Interpreter', 'Latex');

% Set x-axis ticks and labels
xticks(Nvec);
set(gca, 'XTickLabel', Nvec, 'FontName', 'Times', 'FontSize', 26);

% Set y-axis ticks in scientific notation
yticks(10.^(-10:2:-2));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('10^{%d}', x), -10:2:-2, ...
    'UniformOutput', false), 'FontName', 'Times', 'FontSize', 26);

% Set axis limits with small margins
set(gca, 'XLim', [Nvec(1)*0.9, Nvec(end)*1.1]);

% Set y-axis limits based on actual error range
all_errors_vector = cell2mat(all_errors);
valid_errors = all_errors_vector(all_errors_vector > 0);
if ~isempty(valid_errors)
    y_lower = max(min(valid_errors)*0.5, 1e-15);
    y_upper = min(max(valid_errors)*2, 1e10);
    set(gca, 'YLim', [y_lower, y_upper]);
end

% Set figure size and position [left, bottom, width, height] in pixels
set(gcf, 'Position', [100, 100, 700, 500]);

% Customize grid appearance
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';  % Dashed grid lines
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray
ax.GridAlpha = 0.7;      % Semi-transparent
ax.LineWidth = 1.5;      % Axis line width
ax.XMinorGrid = 'off';   % No minor grid lines
ax.YMinorGrid = 'off';

% Final plot formatting
box on;  % Add bounding box
hold off;