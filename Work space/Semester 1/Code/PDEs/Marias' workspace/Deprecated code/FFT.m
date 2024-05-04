% https://www.youtube.com/watch?v=BkA7ncY0b7I

clear all, close all, clc

L = 100; % length domain
N = 1000; % n discretization points
dt = L / N; % Changing dx to dt
t = -L / 2:dt:L / 2 - dt; % time domain

% r2
Beta_2 = 0.83;
Beta_1 = 563.3;
a = (-Beta_2 / 2) / 1i;
b = (-Beta_1 / 2) / 1i;


% discrete wavenumbers
j = (2 * pi / L) * [-N / 2:N / 2 - 1];
j = fftshift(j); % re-order

%% Initial conditions
u0 = zeros(size(t)); % Initialize for time domain
%u0((L / 2 - L / 10) / dt:(L / 2 + L / 10) / dt) = 1;
for i=1:N
    u0(i, 1)=exp(-t(i)^2);  % Gaussian pulse
end


%% Fourier Frequency domain
x = 0:0.1:20; % Replacing x
[x, uhat] = ode45(@(x, uhat) rhspde(x, uhat, j, a,b), x, fft(u0));

u = zeros(length(x), length(t)); % Switch dimensions

for k = 1:length(x) % IFFT to return space domain (replacing t)
    u(k, :) = ifft(uhat(k, :));
end

%% Plot
figure;

%{
%%Waterfall
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
%}

%{ 
%%2D standard
% Plot the real part
subplot(2, 1, 1);
plot(t, real(u), 'b');
xlabel('Time (t)');
ylabel('Real(u)');
title('Real part of u(t)');

% Plot the imaginary part
subplot(2, 1, 2);
plot(t, imag(u), 'r');
xlabel('Time (t)');
ylabel('Imag(u)');
title('Imaginary part of u(t)');

% Adjust the plot appearance
sgtitle('Real and Imaginary Parts of u(t)');

% Plot both real and imaginary parts together
plot(t, real(u), 'b', t, imag(u), 'r');
xlabel('Time (t)');
ylabel('Real and Imaginary Parts');
legend('Real(u)', 'Imag(u)');
title('Real and Imaginary Parts of u(t)');
%}

% Plot the modulus
plot(t, abs(u).^2, 'b');
xlabel('Time (t)');
ylabel('(|u|^2)');
title('|u|^2(t)');
