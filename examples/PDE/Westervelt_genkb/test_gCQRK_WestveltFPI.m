% test the variable Runge-Kutta method for the  wave equations with a
% convolution term:c2d(1-2*kap*u)u_{tt}-cl\Deltau-kc*k\ast\Delta(u_t)=2*kap*(u_t)^2+f(x,t);
clc; clf; close all; clear

% Add quadrature functions path
addpath('../../../quadratures');  

% Simulation parameters
Tf = 1;       % Final time
RK = 2;       % Runge-Kutta order
alp_values = [0.25 0.5 0.75]; % Fractional orders to test
Ord =1 ;      % Theoretical convergence rate

% Wave equation coefficients
kc = 1;       % Convolution coefficient 
r = 1;        % Kernel parameter
c2d = 1;      % u_tt coefficient 
cl = 1;       % Laplace coefficient
kap=0;
% Spatial discretization
J = 400;      % Number of spatial intervals
a = 0; b = 1; % Domain [0,1]
h = (b-a)/J;  % Spatial step size
x1 = ((a+h):h:b-h)'; % Spatial grid points

% Initial conditions and forcing term
f = @(t) sin(pi*x1)*(1+log(t)); % Source term
u0 = 0*sin(pi*x1);              % Initial displacement 
v0 = sin(pi*x1);                % Initial velocity

% Error storage
Nvec = 2.^[2:8]';             % Time step sizes to test
all_errors = cell(length(alp_values),1);

%% Main convergence test loop
for alp_idx = 1:length(alp_values)
    alp = alp_values(alp_idx);
    bet = alp;
    grad = max(1,Ord);
    
    % Fractional kernel function
    Kfun = @(z) kc*((z+r*speye(size(z)))^alp\speye(size(z)));
    
    % Reference solution
    [U] = gCQRK_WestveltFPI(u0,f,Kfun,alp,RK,Nvec(1),Tf,grad,J,h,cl,v0,kap,c2d,r);
    plot(h:h:b-h,U(:,end,end));
    % Error computation
    E = [];
    for k = 1:length(Nvec)
        N = Nvec(k);
        % Refined solution
        U_ref = gCQRK_WestveltFPI(u0,f,Kfun,alp,RK,N*2,Tf,grad,J,h,cl,v0,kap,c2d,r);
        % L2 error computation
        eL2 = sum(abs(U(:,end,:)-U_ref(:,end,1:2:end)).^2)*h;
        E = [E max(sqrt(eL2))];
        U = U_ref;
    end
    
    all_errors{alp_idx} = E;
    
    % Display convergence rates
    rates = log2(E(1:end-1)./E(2:end));
    fprintf('\nResults for alp = %.2f:\n', alp);
    disp(['RK: ', num2str(RK)]);
    disp(['alp: ', num2str(alp)]);
    disp(['bet: ', num2str(bet)]);
    disp('Convergence rates:');
    disp(rates);
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
loglog(Nvec, all_errors{2}(4)*(Nvec(4)./Nvec).^Ord, 'k', 'LineWidth', 2);
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