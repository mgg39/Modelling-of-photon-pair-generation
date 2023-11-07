% https://www.youtube.com/watch?v=BkA7ncY0b7I

clear all, close all, clc

L = 100; % length domain
N = 1000; % n discretization points
dt = L / N; % Changing dx to dt
t = -L / 2:dt:L / 2 - dt; % time domain

% r2
Beta_2 = 0.83;
a = (-Beta_2 / 2) / 1i;

% discrete wavenumbers
kappa = (2 * pi / L) * [-N / 2:N / 2 - 1];
kappa = fftshift(kappa); % re-order

%% Initial conditions
u0 = zeros(size(t)); % Initialize for time domain
u0((L / 2 - L / 10) / dt:(L / 2 + L / 10) / dt) = 1;

%% Fourier Frequency domain
x = 0:0.1:20; % Replacing t
[x, uhat] = ode45(@(x, uhat) rhspde(x, uhat, kappa, a), x, fft(u0));

u = zeros(length(x), length(t)); % Switch dimensions

for k = 1:length(x) % IFFT to return space domain (replacing t)
    u(k, :) = ifft(uhat(k, :));
end

%% Plot
figure;

% Real part
subplot(2, 1, 1);
h_real = waterfall(t, x(1:10:end), real(u(1:10:end, :))); % Switch t and x
set(h_real, 'LineWidth', 5, 'FaceAlpha', 0.5);
colormap(flipud(jet) / 1.5);
set(gca, 'FontSize', 32);
xlabel('Time'); % Switch labels
ylabel('Space');
zlabel('Real Part');
title('Real Part of u');
set(gcf, 'Position', [1500 500 1750 1200]);

% Imaginary part
subplot(2, 1, 2);
h_imag = waterfall(t, x(1:10:end), imag(u(1:10:end, :))); % Switch t and x
set(h_imag, 'LineWidth', 5, 'FaceAlpha', 0.5);
colormap(flipud(jet) / 1.5);
set(gca, 'FontSize', 32);
xlabel('Time'); % Switch labels
ylabel('Space');
zlabel('Imaginary Part');
title('Imaginary Part of u');
