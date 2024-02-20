clear all, close all
clc
%% Constants

N = 1024;
c = 299792458;    %Speed of light
Zmax = 0.03;      %length of waveguide

load("photon_disp.mat"); lscan_photon = lamscan; neff_photon=neff;   %Data for i and s pulses
load("pump_disp.mat"); lscan_pump = lamscan; neff_pump=neff;         %Data for pump pulse
lscan_pump=lscan_pump*10^-6;                                         %Converting from um to m

lambda_min = 1.35e-6; % Minimum and maximum wavelength in meters
lambda_max = 1.65e-6; % Maximum and maximum wavelength in meters

omega_max = 2*pi*c/lambda_min; 
omega_min = 2*pi*c/lambda_max;

n_min = min(neff_photon); 
n_max = max(neff_photon);

%diff array lengths to ensure proper size multiplication
wi = linspace(omega_min, omega_max, 1000); 
ws = linspace(omega_min, omega_max, 1000);
wp_0 = wi + ws;
ni = linspace(n_max, n_min, 1000); 
ns = linspace(n_max, n_min, 1000);

[Ws,Wi] = meshgrid(ws,wi);  %Ws vertical, Wi horizontal
[Ns, Ni] = meshgrid(ns, ni);   %Converting arrays to meshgrids 
Wp = Ws + Wi;  %freq_p = freq_s + freq_i


fitType = @(a, x) a ./ x;
startPoint = [1.5e-6]; % Initial guess for the parameter
params = lsqcurvefit(fitType, startPoint, lscan_pump, neff_pump);  %fitting pump data to above equation
lambda_0 = params(1); %defining lambda_0

Np = lambda_0*Wp./(2*pi*c);    %creating meshgrid of n_eff for pump pulse

%% Betas

Beta_s = Ws.*Ns./c; 
Beta_i = Wi.*Ni./c; 
Beta_p = Wp.*Np./c;   %calculating beta values for all 3 pulses
delta_beta = Beta_p - Beta_s - Beta_i;  %calculating delta_beta

phi = my_sinc(delta_beta,Zmax);    
disp(phi);

%% Plotting - M
% Plotting phi as contour plots
figure;

% Wi vs Ws
subplot(1,2,1);
contour(ws, wi, phi, 20, 'LineWidth', 1.5, 'LineColor', 'k');
xlabel('Frequency (s)');
ylabel('Frequency (i)');
title('I');

% Wi vs Wp0
subplot(1,2,2);
contour(wp_0, wi, phi, 20, 'LineWidth', 1.5, 'LineColor', 'k');
xlabel('Frequency (p)');
ylabel('Frequency (i)');
title('S');

%% Plotting - C
figure;

pcolor(Ws, Wi, phi)
colorbar