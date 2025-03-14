%% chebpde_Heat_1.m -- an executable m-file for solving a partial differential equation
% Automatically created in CHEBGUI by user maye.
% Created on February 12, 2025 at 21:49.

%% Problem description.
% Solving
%   u_t = u'' + exp(u),
% for x in [-1,1] and t in [0,3.544], subject to
%   u = 0 at x = -1
% and
%   u = 0 at x = 1

%% Problem set-up
% Create an interval of the space domain...
dom = [-1 1];
%...and specify a sampling of the time domain:
t = 0:.001:2;

% Make the right-hand side of the PDE.
pdefun = @(t,x,u) diff(u,2)+exp(u);

% Assign boundary conditions.
bc.left = @(t,u) u;
bc.right = @(t,u) u;

% Construct a chebfun of the space variable on the domain,
x = chebfun(@(x) x, dom);
% and of the initial condition.
u0 = chebfun(0,dom);

%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-6, 'Ylim', [0,8]);

%% Call pde15s to solve the problem.
[t, u] = pde15s(pdefun, t, u0, bc, opts);

%% Plot the solution.
%waterfall(u, t)
%xlabel('x'), ylabel('t')

xx = linspace(dom(1), dom(2), 200);
Uvals = zeros(length(xx), length(t));

for k = 1:length(t)
    uk = u(:, k);
    Uvals(:, k) = uk(xx);
end

%% Heat figure
figure;
imagesc(xx, t, Uvals');  
set(gca, 'YDir', 'normal');  
colorbar;
xlabel('x');
ylabel('t');