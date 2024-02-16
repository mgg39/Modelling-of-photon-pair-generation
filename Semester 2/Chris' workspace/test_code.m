clear all, close all
clc
%% Constants

N = 1024;
c = 299792458;
Zmax = 0.03;

load("photon_disp.mat"); lscan_photon = lamscan; neff_photon=neff;
load("pump_disp.mat"); lscan_pump = lamscan; neff_pump=neff;

%% Frequencies 
% Generate array of frequencies
lambda_min = 1.35e-6; lambda_max = 1.65e-6; % Minimum and maximum wavelength in meters
omega_max = 2*pi*c/lambda_min; omega_min = 2*pi*c/lambda_max;

n_min = min(neff_photon); n_max = max(neff_photon);

%diff to ensure proper size multiplication
wi = linspace(omega_min, omega_max, 1000); ws = linspace(omega_min, omega_max, 1020);
ni = linspace(n_max, n_min, 1000); ns = linspace(n_max, n_min, 1020);

[Ws,Wi] = meshgrid(ws,wi); [Ns, Ni] = meshgrid(ns, ni);

%% Betas

Beta_s = Ws.*Ns./c; Beta_i = Wi.*Ni./c; Beta_p = ((Ws+Wi).*4*pi^2*c).^-1;
delta_beta = Beta_p - Beta_s - Beta_i;

phi = sinc(delta_beta, Zmax);

%% Plotting

pcolor(Ws, Wi, phi)
colorbar
