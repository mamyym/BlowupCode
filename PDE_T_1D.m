clear; close all; clc;

% Parameters 
beta = 1;       % beta = 
C_new = 1;       % C_new = , C_new = a*C_old
D_T  = 1;     % D_T = 

% Let T0 vary within a certain range
T0_min = -1;    
T0_max = 3;   % 24.655 
N = 1000;        % number of points
T0_values = linspace(T0_min, T0_max, N);

% Compute Tf based on the formula
Tf_values = T0_values - (2/beta) .* ...
    log( cosh( sqrt(beta*C_new/(2*D_T)) .* exp( (beta .* T0_values)/2 ) ) );

% Plot Bifurtion Diagram
plot(Tf_values, T0_values, 'LineWidth', 2, 'DisplayName', 'Bifurcation Curve');
xlabel('T_f', 'FontSize', 16);
ylabel('T_0', 'FontSize', 16);
set(gca, 'FontSize', 16);
grid on;

hold on;

% Find the maximum T_f (T_c) and corresponding T_0 (T0_c) ---
[T_c, idx] = max(Tf_values);  % T_c is the max value; idx is its index
T0_c = T0_values(idx);        % The T0 that corresponds to T_c

% Mark this point on the plot ---
plot(T_c, T0_c, 'ro', 'MarkerSize', 8, 'LineWidth', 2, ...
    'DisplayName', 'Critical Point');  
text(T_c - 1.0, T0_c, sprintf('T_c=%.2f, T_0=%.2f \', T_c, T0_c), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 16);

% Draw a vertical line at x = T_c
xline(T_c, '--r', 'LineWidth', 2, 'DisplayName', 'Critical Line');

legend('Location','best', 'FontSize', 12);

hold off;

% open another new figure
figure; 
plot(T0_values, Tf_values, 'LineWidth', 2, 'DisplayName', 'Bifurcation Curve');
hold on; grid on;

xlabel('T_0', 'FontSize', 14);
ylabel('T_f', 'FontSize', 14);
set(gca, 'FontSize', 12);

plot(T0_c, T_c, 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Critical Point');

text(T0_c, T_c + 0.5, sprintf('T_c=%.2f, T_0=%.2f', T_c, T0_c), ...
    'HorizontalAlignment', 'center', ...   
    'VerticalAlignment', 'bottom', ...     
    'FontSize', 16);

yline_obj = yline(T_c, '--r', ...
    'LabelHorizontalAlignment','left', ...
    'LabelVerticalAlignment','top', ...
    'DisplayName','Critical Line');
set(yline_obj, 'LineWidth', 1.5);

legend('Location','best', 'FontSize', 12);
hold off;
