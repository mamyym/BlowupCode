%% 
beta = 1;     
C_new = 1;      
D_T_values = linspace(0.01, 10, 10000);

T0_min = -20;
T0_max = 20;
N = 10000;
T0_values = linspace(T0_min, T0_max, N);

Tc_values = zeros(size(D_T_values));

%% 
for i = 1:length(D_T_values)
    D_T = D_T_values(i);

    Tf = T0_values - (2/beta) * log( cosh( sqrt(beta * C_new/(2*D_T)) .* exp((beta * T0_values)/2) ) );

    Tc_values(i) = max(Tf);
end

%% 
figure;
plot(D_T_values, Tc_values, 'b-', 'LineWidth', 2);
xlabel('D_T', 'FontSize', 16);
ylabel('T_c', 'FontSize', 16);
set(gca, 'FontSize', 16);
legend('T_c', 'Location', 'best', 'FontSize', 14);
grid on;
