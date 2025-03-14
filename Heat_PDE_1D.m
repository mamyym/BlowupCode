%clear; close all; clc;

% Spatial domain [xL, xR] = [-1, 1]
xL = -1.0;
xR =  1.0;
Nx = 401;                  % Number of spatial grid points (including boundaries)
x  = linspace(xL, xR, Nx);
dx = x(2) - x(1);         % Spatial step size

% Time domain [0, tmax]
tmax = 1.0;              
dt   = 1e-5;          % Explicit scheme requires a sufficiently small time step
Nt   = round(tmax/dt); % Number of time steps
t    = linspace(0, tmax, Nt+1);

% PDE coefficients
D_T    = 1.0;  
c_new  = 1.0;
beta   = 1.0;

% Boundary and initial conditions
Tf     = 0.0;  % Boundary temperature: T(-1,t) = T(1,t) = Tf
T0     = 0.0;  % Initial temperature: T(x,0) = T0

%% 2. Allocate arrays and set initial/boundary values
% T_data(i, n) represents the temperature at spatial grid point i and time step n
T_data = zeros(Nx, Nt+1);

% Initial condition: Temperature distribution in space at t = 0
T_data(:, 1) = T0;

% Boundary conditions: x = -1 corresponds to i = 1, x = 1 corresponds to i = Nx
T_data(1,   :) = Tf;  
T_data(Nx,  :) = Tf;

%% 3. Explicit finite difference main loop
for n = 1:Nt
    for i = 2:(Nx-1)
        % Second-order spatial derivative (Txx)
        Txx = (T_data(i+1,n) - 2*T_data(i,n) + T_data(i-1,n)) / dx^2;
        % Explicit scheme update
        T_data(i, n+1) = T_data(i,n) ...
                         + dt * ( D_T * Txx ...
                                  + c_new * exp(beta * T_data(i,n)) );
    end
    % Keep boundary temperature unchanged
    T_data(1,   n+1) = Tf;
    T_data(Nx,  n+1) = Tf;
end

%% 4. Plot the heatmap of T(x,t)
figure;
imagesc(x, t, T_data');    % Transpose so that rows correspond to t, columns to x
set(gca, 'YDir', 'normal');
colorbar;
colormap('hot');
xlabel('x');
ylabel('t');