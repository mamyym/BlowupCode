%% 
betas   = [0.5, 1, 2, 5, 10];   
C_new   = 1;                 
D_T     = 5;                 
T0_min  = -30;                 
T0_max  = 24;              
N       = 1000;          
T0_values = linspace(T0_min, T0_max, N);

%% 
figure;
hold on; grid on;
colors = lines(length(betas)); 

for i = 1:length(betas)
    beta = betas(i);
    
    Tf_values = T0_values - (2/beta) .* log( cosh( sqrt(beta*C_new/(2*D_T)) ...
                   .* exp((beta .* T0_values)/2) ) );

    [T_c, idx] = max(Tf_values);
    T0_c       = T0_values(idx);

    plot(T0_values, Tf_values, 'LineWidth', 2, 'Color', colors(i,:), ...
         'DisplayName', sprintf('\\beta = %.2g', beta));

    plot(T0_c, T_c, 'o', 'MarkerSize', 5, 'LineWidth', 2, ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
         'HandleVisibility','off');

    yline(T_c, '--', 'Color', colors(i,:), 'LineWidth', 1.5, ...
          'HandleVisibility','off');
    
end

xlabel('T_0', 'FontSize', 16);
ylabel('T_f', 'FontSize', 16);
ylim([-10, 10]);
set(gca, 'FontSize', 16);
legend('Location','best', 'FontSize', 16);

hold off;
