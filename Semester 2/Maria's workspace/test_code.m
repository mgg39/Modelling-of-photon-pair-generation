clear all, close all
clc

% Constants
N = 1024;
c = 299792458; %Speed of light
Zmax = 0.03; %length of waveguide

load("photon_disp.mat"); 
lscan_photon = lamscan; 
neff_photon=neff;   %Data for i and s pulses
lscan_photon=lscan_photon*10^-6;   %Converting from um to m

load("pump_disp.mat"); lscan_pump = lamscan; neff_pump=neff; %Data for pump pulse
lscan_pump=lscan_pump*10^-6;   %Converting from um to m

lambda_min = 1.35e-6; % Minimum and maximum wavelength in meters
lambda_max = 1.65e-6; % Maximum and maximum wavelength in meters

omega_max = 2*pi*c/lambda_min; 
omega_min = 2*pi*c/lambda_max;

wi = linspace(omega_min, omega_max, 1000); 
ws = linspace(omega_min, omega_max, 1000);
wp_0 = wi + ws;

wscan_photon=2*pi*c./(lscan_photon);

ns=spline(wscan_photon,neff_photon,ws);
ni=spline(wscan_photon,neff_photon,wi);

wscan_pump=2*pi*c./(lscan_pump); 

[Ws,Wi] = meshgrid(ws,wi);  %Ws vertical, Wi horizontal
[Ns, Ni] = meshgrid(ns, ni);   %Converting arrays to meshgrids 
Wp = Ws + Wi;  %freq_p = freq_s + freq_i

Np = spline(wscan_pump,neff_pump,Ws+Wi);

% Betas
Beta_s = Ws.*Ns./c; 
Beta_i = Wi.*Ni./c; 
Beta_p = Wp.*Np./c;   %calculating beta values for all 3 pulses

delta_beta = Beta_p - Beta_s - Beta_i;  %calculating delta_beta
phi = my_sinc(delta_beta.*Zmax/2);

w0=2*pi*c/750e-9;
alpha=exp(-(Ws+Wi-w0).^2*0.25e-24/2); %TODO: change - assumed gaussian not sech

figure
pcolor(Ws,Wi,alpha);
shading interp;
xlabel('Frequency (p)');
ylabel('Frequency (i)');
title('Alpha');

figure
pcolor(Ws,Wi,alpha.*phi);
shading interp;
xlabel('Frequency (p)');
ylabel('Frequency (i)');
title('Phi');

%---------------------------------------------------------------
% Integral
%---------------------------------------------------------------

%placeholder vars
p = 1; %gamma from ODEs

% placeholder grid domain vals
z = 0:.1:3; 

%currently broken due to array size differences ----------------
disp(size(Beta_p));
disp(size(z));
%---------------------------------------------------------------

f = p * exp(1i * Beta_p * z);

I = trapz(f);
disp(I);