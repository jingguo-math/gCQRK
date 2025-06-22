% Test the Runge-Kutta based gCQ method for the fractional integral
clc
clf
close all
clear all
% Add path to quadrature functions (relative path from current location)
addpath('../../quadratures');
%% Numerical Parameters
RK = 2;                 % Runge-Kutta method order (2nd order method)
Ord = 3;                % Expected convergence order
Tf = 1;                 % Final simulation time

%% Problem Parameters
alp = 0.5;              % Fractional order of integration
bet_values = [-0.1 0.3 0.7 1.1 1.5]; % Test different regularity exponents bet
Nvec = 2.^(3:10);    % Vector of grid sizes (N = [8, 16, 32, ..., 1024])

% Kernel function for fractional integral: K(z) = z^(-alp)
Kfun = @(z) z.^(-alp);

%% Initialize storage for errors
all_errors = cell(length(bet_values), 1); % Cell array to store errors for each 

%% Main Convergence Study Loop
for bet_idx = 1:length(bet_values)
    bet = bet_values(bet_idx);
    
    %% Compute Optimal Grading Parameter
    % The grading parameter ensures proper resolution near t=0
    % Formula: grad = max(1, p/(alp+bet)) where p is the scheme order
    grad = max(1, Ord/(alp + bet));
    
    %% Define Exact Solution and Right-Hand Side
    % Right-hand side function: f(t) = t^bet
    f = @(t) t.^bet;
    
    % Exact solution of fractional integral
    sol = @(t) gamma(bet+1) * t.^(alp+bet) / gamma(alp+bet+1);
    
    %% Error Computation for Different Grid Sizes
    E = []; % Initialize error vector
    
    for j = 1:length(Nvec)
        N = Nvec(j);
        
        % Call the gCQ Runge-Kutta solver
        [C, eamax] = cqrk_varn0(sol, f, Kfun, RK, N, Tf, alp, grad);
        
        % Store maximum absolute error
        E = [E, eamax];
    end
    
    % Store errors for this bet value
    all_errors{bet_idx} = E;
    
    %% Compute and Display Convergence Rates
    rate = log2(E(1:end-1)./E(2:end)); % log2 of error ratios
    
    disp('RK method:');
    disp(RK);
    disp('alp:');
    disp(alp);
    disp('bet:');
    disp(bet);
    disp('Convergence rates:');
    disp(rate);
end

%% Visualization of Results
figure(1);

% Plot styling parameters
line_styles = {'d-', 'p:', '*--', 's-.', 'o--'}; % Diamond, pentagon, asterisk, square, circle
colors = [
    0 0.4470 0.7410;    % Blue
    0.8500 0.3250 0.0980; % Orange
    0.9290 0.6940 0.1250; % Yellow
    0.4940 0.1840 0.5560; % Purple
    0.4660 0.6740 0.1880  % Green
    ];

%% Plot Reference Convergence Line
% Shows ideal slope of -Ord (expected convergence rate)
loglog(Nvec, all_errors{2}(3)*(Nvec(3)./Nvec).^Ord, 'k', 'LineWidth', 2);
hold on

%% Plot Error Curves for Each bet Value
legend_entries = cell(1, length(bet_values) + 1);
legend_entries{1} = sprintf('$slope=%.0f$', -Ord); % Reference line legend

for bet_idx = 1:length(bet_values)
    % Create legend entry for this bet value
    legend_entries{bet_idx+1} = sprintf('$\\beta=%.1f$', bet_values(bet_idx));
    
    % Plot error curve with specified style
    loglog(Nvec, all_errors{bet_idx}, line_styles{bet_idx}, ...
        'LineWidth', 2, 'MarkerSize', 12, 'Color', colors(bet_idx,:));
    hold on
end

%% Figure Formatting
% Axis labels with LaTeX formatting
xlabel('$N$ (Number of time steps)', 'FontSize', 30, 'Interpreter', 'Latex');
ylabel('Maximum Absolute Error', 'FontSize', 30, 'Interpreter', 'Latex');

% Legend configuration
legend(legend_entries, 'Location', 'southwest', 'FontSize', 23, 'Interpreter', 'Latex');

% X-axis ticks
xticks(Nvec);
set(gca, 'XTickLabel', Nvec, 'FontName', 'Times', 'FontSize', 26);

% Y-axis ticks in scientific notation
yticks(10.^(-10:2:-2));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('10^{%d}', x), -10:2:-2, ...
    'UniformOutput', false), 'FontName', 'Times', 'FontSize', 26);

% Axis limits
set(gca, 'XLim', [Nvec(1)*0.9, Nvec(end)*1.1]); % X-axis with 10% margin

% Dynamic y-axis limits based on error range
all_errors_vector = cell2mat(all_errors);
valid_errors = all_errors_vector(all_errors_vector > 0);
if ~isempty(valid_errors)
    y_lower = max(min(valid_errors)*0.5, 1e-15);
    y_upper = min(max(valid_errors)*2, 1e10);
    set(gca, 'YLim', [y_lower, y_upper]);
end

% Figure size and position [left, bottom, width, height] in pixels
set(gcf, 'Position', [100, 100, 700, 500]);

% Grid customization
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';      % Dashed grid lines
ax.GridColor = [0.7, 0.7, 0.7]; % Light gray
ax.GridAlpha = 0.7;           % Semi-transparent
ax.LineWidth = 1.5;           % Axis line width
ax.XMinorGrid = 'off';        % No minor x-grid
ax.YMinorGrid = 'off';        % No minor y-grid

box on;  % Add bounding box
hold off;