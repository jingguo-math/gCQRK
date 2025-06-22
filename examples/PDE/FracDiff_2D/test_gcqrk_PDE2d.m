%% gCQRK for 2d PDE D_t^\alp u-c_u*Lap*u=f(f is a given function) at prescribed general time mesh tn, n=1,...,N
clc
clf
close all
clear

% Add path to quadrature functions (three levels up from current location)
addpath('../../../quadratures');  

% Numerical method parameters
RK=2;       % Runge-Kutta method order (2nd order)
Ord=3;      % Expected convergence order
Tf=1;       % Final simulation time

% Problem parameters
alp_values = [0.25 0.5 0.75]; % Fractional derivative orders to test (only 0.75 in this case)

% Spatial domain setup
lps=1;      % Half domain length [-lps, lps]
a=-lps;     % Left boundary
b=lps;      % Right boundary
J=256;       % Number of spatial intervals
h=(b-a)/J;  % Spatial step size
E=[];       % Initialize error array
c_u=1;      % Coefficient for Laplace operator

% Spatial grid points (excluding boundaries)
x1=((a+h):h:b-h)';  % x-coordinates
y1=x1';             % y-coordinates (same as x for square domain)

% Time discretization parameters
Nvec=2.^[4:10]';   
maxIt=length(Nvec); % Number of convergence tests
E= zeros(maxIt,1);  % Initialize error storage

% Initialize cell array to store errors for each alpha value
all_errors = cell(length(alp_values), 1);

%% Main loop over fractional orders (α values)
for alp_idx = 1:length(alp_values)
    alp = alp_values(alp_idx);
    bet=alp;        % Solution regularity exponent matches fractional order
    grad=max(1,Ord/bet); % Grading parameter for time mesh
    
    E=[]; % Initialize error vector for current alpha
    
    % Exact solution: u(x,y,t) = cos(πx/2)cos(πy/2) * t^β
    sol = @(t) cos(pi/2*x1)*cos(pi/2*y1).*t.^bet;
    
    % Right-hand side function f(x,y,t):
    % Contains fractional derivative term and Laplace operator term
    f=@(t) cos(pi/2*x1)*cos(pi/2*y1)*gamma(bet+1).*t.^(bet-alp)/gamma(bet-alp+1)+...
        pi^2/2*cos(pi/2*x1)*cos(pi/2*y1).*c_u.*t.^bet;
    
    %% Convergence study over different time step counts
    for k = 1:maxIt
        % Call 2D fractional PDE solver:
        % - sol: exact solution
        % - f: right-hand side
        % - RK: Runge-Kutta order
        % - Nvec(k): current time step count
        % - Tf: final time
        % - grad: time mesh grading
        % - J: spatial intervals
        % - alp: fractional order
        % - h: spatial step size
        % - c_u: Laplace coefficient
        [C,eL2max] = cqrk_varPDE2d(sol,f,RK,Nvec(k),Tf,grad,J,alp,h,c_u);
        
        E(k,1)=eL2max; % Store maximum L2 error
    end
    
    % Store errors for this alpha value
    all_errors{alp_idx} = E;
    
    % Compute convergence rates (log2 of error ratios between successive N)
    rate = log2(E(1:end-1)./E(2:end));
    
    %% Display convergence results
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
line_styles = {'d-', 'p:', '*--', 's-.', 'o--'}; % Diamond, pentagon, asterisk, square, circle
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

% Create legend entries (reference line + alpha values)
legend_entries = cell(1, length(alp_values) + 1);
legend_entries{1} = sprintf('$slope=%.0f$', -Ord); % Reference slope entry

% Plot error curves for each alpha value
for alp_idx = 1:length(alp_values)
    legend_entries{alp_idx+1} = sprintf('$\\alpha=%.2f$', alp_values(alp_idx));
    loglog(Nvec, all_errors{alp_idx}, line_styles{alp_idx}, ...
        'LineWidth', 2, 'MarkerSize', 12, 'Color', colors(alp_idx,:));
    hold on
end

%% Figure formatting
% Axis labels with LaTeX formatting
xlabel('$N$', 'FontSize', 30, 'Interpreter', 'Latex');
ylabel('Maximum Discrete $L^2$ Norm Error', 'FontSize', 30, 'Interpreter', 'Latex');

% Legend configuration
legend(legend_entries, 'Location', 'southwest', 'FontSize', 23, 'Interpreter', 'Latex');

% X-axis ticks and labels
xticks(Nvec);
set(gca, 'XTickLabel', Nvec, 'FontName', 'Times', 'FontSize', 26);

% Y-axis ticks in scientific notation
yticks(10.^(-10:2:-2));
set(gca, 'YTickLabel', arrayfun(@(x) sprintf('10^{%d}', x), -10:2:-2, ...
    'UniformOutput', false), 'FontName', 'Times', 'FontSize', 26);

% Axis limits with small margins
set(gca, 'XLim', [Nvec(1)*0.9, Nvec(end)*1.1]);

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