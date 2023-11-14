clear all, close all, clc

L = 100; % length domain
N = 1000; % n discretization points
dt = L / N; % Changing dx to dt
t = -L / 2:dt:L / 2 - dt; % time domain

% r2
Beta_2 = (0.83)*(10.^(-24));
Beta_1 = (563.3)*(10.^(-12));
a = (-Beta_2 / 2) / 1i;
b = (-Beta_1 / 2) / 1i;

% discrete wavenumbers
j = (2 * pi / L) * [-N / 2:N / 2 - 1];
j = fftshift(j); % re-order

%% Initial conditions
u0 = zeros(size(t)); % Initialize for time domain
%u0((L / 2 - L / 10) / dt:(L / 2 + L / 10) / dt) = 1;
for i=1:N
    u0(i) = sech(t(i));  % Gaussian pulse
end

%% Fourier Frequency domain
x = 0:0.1:20; % Replacing x
[x, uhat] = ode45(@(x, uhat) rhspde(x, uhat, j, a, b), x, fft(u0));

u = zeros(length(x), length(t)); % Switch dimensions

for k = 1:length(x) % IFFT to return space domain (replacing t)
    u(k, :) = ifft(uhat(k, :));
end


%% Plot
figure;

pcolor(t,x,abs(u).^2)
shading interp
xlabel('t(ns)')
ylabel('x(m)')