clear all, close all
clc
%% Constants

N = 1024;

load("photon_disp.mat");
lscan_photon = lamscan; neff_photon=neff;
load("pump_disp.mat");
lscan_pump = lamscan; neff_pump=neff;

%% Frequencies 
% Generate array of frequencies
lambda_min = 1.35e-9; % Minimum wavelength in meters
lambda_max = 1.65e-9; % Maximum wavelength in meters

%diff to ensure proper size multiplication
wi = wavelength_to_frequency_wi(lambda_min, lambda_max);
ws = wavelength_to_frequency_ws(lambda_min, lambda_max);

[S,I] = meshgrid(wi,ws);

wi_matrix = meshgrid(wi,wi);

num_columns = length(ws); 

[X, ~] = meshgrid(1:num_columns, 1:numel(ws));

ws_matrix = repmat(ws, 1, num_columns);

disp(ws_matrix);


%% Betas
%Betas = interp1(x, lscan_photon, xq, linear);


%% Plotting
%pcolor(S, I, phi)
