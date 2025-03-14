%% 
beta = 1;    
D_T  = 1;    
C_new_values = linspace(0.01, 10, 10000);

T0_min = -20;
T0_max = 20;
N = 1000;
T0_values = linspace(T0_min, T0_max, N);

Tc_values = zeros(size(C_new_values));

%% 
for i = 1:length(C_new_values)
    C_new = C_new_values(i);

    Tf = T0_values - (2/beta) * log( cosh( sqrt(beta * C_new/(2*D_T)) * exp((beta * T0_values)/2) ) );

    Tc_values(i) = max(Tf);
end

%% 
figure;
plot(C_new_values, Tc_values, 'b-', 'LineWidth', 2);
xlabel('c_{new}', 'FontSize', 16);
ylabel('T_c', 'FontSize', 16);
set(gca, 'FontSize', 16);
legend('T_c', 'Location', 'best', 'FontSize', 14);
grid on;
