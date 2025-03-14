%% 
C_new_values = [0.1, 0.5, 1, 2, 5];  
beta         = 1;                  
D_T          = 1;                   

T0_min  = -20;                      
T0_max  = 20;                       
N       = 1000;                    
T0_values = linspace(T0_min, T0_max, N);

%% 
figure;
hold on; grid on;
colors = lines(length(C_new_values));  

for i = 1:length(C_new_values)
    C_new = C_new_values(i);

    Tf_values = T0_values - (2/beta) .* log( cosh( sqrt(beta * C_new/(2*D_T)) .* exp((beta .* T0_values)/2) ) );

    [T_c, idx] = max(Tf_values);
    T0_c       = T0_values(idx);

    plot(T0_values, Tf_values, 'LineWidth', 2, 'Color', colors(i,:), ...
         'DisplayName', sprintf('c_{new} = %.2g', C_new));

    plot(T0_c, T_c, 'o', 'MarkerSize', 4, 'LineWidth', 2, ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'HandleVisibility','off');

    yline(T_c, '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility','off');
    
end

xlabel('T_0', 'FontSize', 16);
ylabel('T_f', 'FontSize', 16);
set(gca, 'FontSize', 16);
ylim([-10, 10]);

legend('Location','best', 'FontSize', 14);
hold off;
