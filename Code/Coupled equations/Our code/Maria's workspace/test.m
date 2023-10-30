clear all
clc

% Numerical values for parameters and initial conditions
r2 = 1.0;
s1 = 2.0;
s2 = 3.0;
g = 4.0;
p1 = 5.0;
p2 = 6.0;
F0 = 0.0;
S0 = 0.0;
P0 = 0.0;

% Constants
k = physconst('Boltzmann');

% Symbolic variables
syms F(z,t) S(z,t) P(z,t)

% Convert symbolic equations into numerical functions
UF = @(z,t,F,dFdt) i*(dFdt) - r2*(dFdt - dFdt(2:end-1)) - S.*conj(F).*exp(1i*k*z);
US = @(z,t,S,dSdt) i*(dSdt) - s1*(dSdt) - s2*(dSdt - dSdt(2:end-1)) - (F.^2)/2.*exp(-1i*k*z) + g*P;
UP = @(z,t,P,dPdt) i*(dPdt) - p1*(dPdt) - p2*(dPdt - dPdt(2:end-1)) + g*S;

% Initial conditions
tspan = [0 10];  % Time span for integration
zspan = [0 1];   % Spatial domain

F0 = @(z) F0;  % You can define your own initial conditions
S0 = @(z) S0;
P0 = @(z) P0;

% Solving the differential equations using ode45
[t, y] = ode45(@(t, y) [UF(z, t, y(1), y(2)); US(z, t, y(1), y(2), y(3)); UP(z, t, y(1), y(2), y(3))], tspan, [F0(z); S0(z); P0(z)], odeset('RelTol', 1e-5));

% Plot 
figure
plot(t,abs(y).^2)