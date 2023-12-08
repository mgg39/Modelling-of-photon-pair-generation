clear all;
close all;
clc;

T = 20; % time domain width
N = 1024; % N discretization points
dt = T / N;
t = linspace(-T/4, 3*T/4 - dt, N)'; % time domain in ps (ps determined by constants)

delta = (2 * pi / T) * linspace(-1, 1, N);

% Constants
Beta_f2 = 0.83e-2; % Units in ps^2/cm
Beta_s2 = 0.22e-2; % Units in ps^2/cm
Beta_p2 = 2.53e-2; % Units in ps^2/cm
Beta_f1 = 563.3e-2; % Units in ps/cm
Beta_s1 = 533.3e-2; % Units in ps/cm
kappa =  6.9e3; % Units in 1/cm
gamma = 1 * 10^1.5; % Units in 1/(cm*sqrt(kW))

C = 1; % Units in 1/cm, C=2 corresponds to a rail separation of x=200nm  

%% Initial conditions

u0 = zeros(3 * N, 1); % Defining an array to represent all pulses together
                       % F pulse represented by first N points, S represented by
                       % N+1 to 2N point, P represented by 2N+1 to 3N points

pulsewidth = 1; % Pulse width of laser (timeframe already scaled to ps with constants)
A = 1; % Amplitude of laser pulse in kiloWatts (kW scaled by constants)
ratio = 2 * asech(1/2) / pulsewidth; % Finding the ratio between the desired pulsewidth and FWHM of a sech curve to scale t by

for c = 1:N
    u0(c) = sqrt(A) * sech(t(c) * ratio) * (1 + 1i) / sqrt(2); % Sech represents laser pulses very well
    % Actual height of sech curve is A^2 (due to abs(u3).^2) so
    % it is scaled by sqrt(A)
end 

% Other pulses remain at 0 for initial conditions

opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-4);

%% Damage threshold

D_T = 14e-3;  % Damage threshold of LiNbO3, units of kJ/cm^2
w = 664e-7;    % Width of waveguide in cm
h = 330e-7;    % Height of waveguide in cm
theta = 70;    % Slant of waveguide wall in degrees

face = w * h + h^2 / tan(theta); % Area of waveguide face in cm^2

% Integral of a*sech(b*x) from -infinity to infinity = a*pi/b
E = sqrt(A) * pi / (ratio * 10^12);

max_PW = 2 * asech(1/2) * D_T * face / (pi * sqrt(A) * 10^-12); % Maximum allowed pulse width for selected pulse strength without damaging crystal
max_PA = (D_T * face * ratio * 10^12 / pi)^2;  % Maximum allowed pulse strength for selected pulse width without damaging crystal

if (E / face) > D_T
    fprintf("The input laser has exceeded the damage threshold of Lithium Niobate \n");
    fprintf("Try changing the pulsewidth to be less than %.3f ps or the pulse strength to be less than %.3f kW \n", max_PW, max_PA);
    return;
end

%% Fourier Frequency domain
zend = 5 * T / (8 * Beta_f1);
z = linspace(0, zend, N)'; % Spatial domain in cm (cm determined by constants)

[z, uhat] = ode45(@(z, uhat) CoupledPDEs(z, uhat, N, delta, Beta_f1, Beta_f2, Beta_s1, Beta_s2, Beta_p2, kappa, gamma, C), z, u0, opts);

u1 = uhat(:, 1:N); % Breaking apart final matrix into 3 respective pulses
u2 = uhat(:, N+1:2*N);
u3 = uhat(:, 2*N+1:3*N);

%% Plot
figure;

subplot(1, 3, 1) % Plotting F pulse
pcolor(t, z, abs(u1).^2)
shading interp
hold on
xlabel('t (ps)')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Pulse power (kW)", "fontsize", 10, "rotation", 270)
title("F")
set(gca, 'TickDir', 'out'); 

subplot(1, 3, 2) % Plotting S pulse
pcolor(t, z, abs(u2).^2)
shading interp
hold on
xlabel('t (ps)')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Pulse power (kW)", "fontsize", 10, "rotation", 270)
title("S")
set(gca, 'TickDir', 'out'); 

subplot(1, 3, 3) % Plotting P pulse
pcolor(t, z, abs(u3).^2)
shading interp
hold on
xlabel('t (ps)')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Pulse power (kW)", "fontsize", 10, "rotation", 270)
title("P")
set(gca, 'TickDir', 'out'); 

%% Peak finder

[Pmax, Idx] = max(abs(u3(:)).^2);
[PmaxRow, PmaxCol] = ind2sub(size(abs(u3).^2), Idx);
Zmax = z(PmaxRow);
Tmax = t(PmaxCol);

fprintf("For an input laser of power %.2f kW and pulsewidth %.1f ps, the Pump pulse has a maximum amplitude of %.2f kW at z = %.2f cm and t = %.1f ps\n", A, pulsewidth, Pmax, Zmax, Tmax)

%% Fourier Transform at each point z

U1_hat = fft(u1, N, 2); % Fourier transform of F pulse
U2_hat = fft(u2, N, 2); % Fourier transform of S pulse
U3_hat = fft(u3, N, 2); % Fourier transform of P pulse

%% Plot Fourier Transforms

figure;

subplot(1, 3, 1) % Plotting Fourier transform of F pulse
pcolor(delta, z, abs(fftshift(U1_hat, 2)).^2)
shading interp
xlabel('Frequency (\delta)')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Field Spectrum", "fontsize", 10, "rotation", 270)
title("F")
set(gca, 'TickDir', 'out');

subplot(1, 3, 2) % Plotting Fourier transform of S pulse
pcolor(delta, z, abs(fftshift(U2_hat, 2)).^2)
shading interp
xlabel('Frequency (\delta)')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Field Spectrum", "fontsize", 10, "rotation", 270)
title("S")
set(gca, 'TickDir', 'out');

subplot(1, 3, 3) % Plotting Fourier transform of P pulse
pcolor(delta, z, abs(fftshift(U3_hat, 2)).^2)
shading interp
xlabel('Frequency (\delta)')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Field Spectrum", "fontsize", 10, "rotation", 270)
title("P")
set(gca, 'TickDir', 'out');

