clear all, close all
clc
%% Constants

N = 1024;
c = 299792458;    %Speed of light
Zmax = 0.03;      %length of waveguide

load("photon_disp.mat"); lscan_photon = lamscan; neff_photon=neff;   %Data for i and s pulses
load("pump_disp.mat"); lscan_pump = lamscan; neff_pump=neff;         %Data for pump pulse
lscan_pump=lscan_pump*10^-6;                                         %Converting from um to m

%% Frequencies 
% Generate array of frequencies, wavelegnths and n_effectives
lambda_min = 1.35e-6; lambda_max = 1.65e-6; % Minimum and maximum wavelength in meters
omega_max = 2*pi*c/lambda_min; omega_min = 2*pi*c/lambda_max;

n_min = min(neff_photon); n_max = max(neff_photon);

%diff array lengths to ensure proper size multiplication
wi = linspace(omega_min, omega_max, 1000); ws = linspace(omega_min, omega_max, 1020);
ni = linspace(n_max, n_min, 1000); ns = linspace(n_max, n_min, 1020);

[Ws,Wi] = meshgrid(ws,wi); [Ns, Ni] = meshgrid(ns, ni);   %Converting arrays to meshgrids 
Wp = Ws + Wi;                                             %freq_p = freq_s + freq_i

fitType = fittype('a/x');                                             %n=lambda_0/lambda
params = fit(lscan_pump, neff_pump, fitType, StartPoint = [1.5e-6]);  %fitting pump data to above equation
lambda_0 = params(1);                                                 %defining lambda_0

Np = lambda_0*Wp./(2*pi*c);    %creating meshgrid of n_eff for pump pulse

%% Betas

Beta_s = Ws.*Ns./c; Beta_i = Wi.*Ni./c; Beta_p = Wp.*Np./c;   %calculating beta values for all 3 pulses
delta_beta = Beta_p - Beta_s - Beta_i;                        %calculating delta_beta

phi = sinc(delta_beta, Zmax);     %phi = sinc(delta_beta*Z/2)

%% Plotting

pcolor(Ws, Wi, phi)
colorbar
