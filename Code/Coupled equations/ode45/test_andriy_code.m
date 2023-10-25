% Clear the workspace and command window
clear all
clc

% Set values for gamma and kappa
gamma = 1;
kappa = 0.1;

% Initial condition for the differential equation
xini = [1, 0];

% Configure options for the ODE solver
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Set the maximum value for the independent variable z
zmax = 50;

% Use the ODE solver (ode45) to solve the differential equation defined by 'rhs'
[z, x] = ode45(@(z, x) rhs(z, x, gamma, kappa), [0, zmax], xini, options);

% Create a new figure for plotting
figure

% Plot the absolute value squared of the solution 'x' against 'z'
plot(z, abs(x).^2)

% Define the right-hand side function for the differential equation
function y = rhs(z, x, gamma, kappa)

% Initialize the output vector
y = zeros(size(x));

% Define the differential equations
y(1) = 1i * gamma * conj(x(1)) * x(2) * exp(1i * kappa * z);
y(2) = 0.5i * gamma * x(1)^2 * exp(-1i * kappa * z);

end

% Add two more variables and initial conditions
yini = [0, 1];

% Solve the extended system of ODEs
[z, y] = ode45(@(z, y) rhs(z, y, gamma, kappa), [0, zmax], [xini, yini], options);

% Create a new figure for plotting
figure

% Plot the solutions for x and y
plot(z, abs(y(:, 1)).^2, z, abs(y(:, 3)).^2)

% Define the right-hand side function for the extended ODE system
function dydz = rhs(z, y, gamma, kappa)

% Extract the variables
x1 = y(1);
x2 = y(2);
y1 = y(3);
y2 = y(4);

% Initialize the output vector
dydz = zeros(4, 1);

% Define the ODEs for the extended system
dydz(1) = 1i * gamma * conj(x1) * x2 * exp(1i * kappa * z);
dydz(2) = 0.5i * gamma * x1^2 * exp(-1i * kappa * z);
dydz(3) = -0.5i * gamma * y1^2 * exp(-1i * kappa * z);
dydz(4) = -1i * gamma * conj(x2) * y1 * exp(1i * kappa * z);
end
