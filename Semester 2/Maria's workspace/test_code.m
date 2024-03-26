clear all, close all
clc

tic;  %Measuring time elapsed over code

%% Constants

N = 2048;
c = 299792458; %Speed of light

load("photon_disp.mat"); %Data for i and s pulses
lscan_photon = lamscan; neff_photon=neff;   
lscan_photon=lscan_photon*10^-6;   %Converting from um to m

load("pump_disp.mat"); lscan_pump = lamscan; neff_pump=neff; %Data for pump pulse 
lscan_pump=lscan_pump*10^-6;   %Converting from um to m

load("PumpPulse.mat");
Zmax = max(z)/100;  %Defining Zmax as the length of the waveguide (converting from cm to m)

%% 

wi = linspace(1.34e15, 1.37e15, 2000); 
ws = linspace(1.146e15, 1.166e15, 2000);

wscan_photon=2*pi*c./(lscan_photon);

ns=spline(wscan_photon,neff_photon,ws);  %using spline to get range of neff that follow same relation between neff_photon and wscan_photon
ni=spline(wscan_photon,neff_photon,wi);  

wscan_pump=2*pi*c./(lscan_pump); 

[Ws,Wi] = meshgrid(ws,wi);     %Converting arrays to meshgrids
[Ns, Ni] = meshgrid(ns, ni);   %Ws vertical, Wi horizontal
Wp = Ws + Wi;                  %freq_p = freq_s + freq_i

Np = spline(wscan_pump,neff_pump,Wp);   %using spline to get range of neff that follow same relation between neff_pump and wscan_pump

%% Betas

Beta_s = Ws.*Ns./c;   %beta = n*w/c
Beta_i = Wi.*Ni./c; 
Beta_p = Wp.*Np./c;   %calculating beta values for all 3 pulses

delta_beta = Beta_p - Beta_s - Beta_i;  %delta_beta = beta_p - beta_s - beta_i
phi = sinc(delta_beta, Zmax);

w0=2*pi*c/750e-9;
alpha=exp(-(Ws+Wi-w0).^2*0.25e-24/2);

%% Plotting

figure
pcolor(Ws,Wi,alpha);
shading interp;
xlabel('\omega_s (Hz)');
ylabel('\omega_i (Hz)');
title('\alpha');
colorbar;
set(gca,'TickDir','out'); 

figure
pcolor(Ws,Wi,phi);
shading interp;
xlabel('\omega_s (Hz)');
ylabel('\omega_i (Hz)');
title('\phi');
colorbar;
set(gca,'TickDir','out'); 

figure
pcolor(Ws,Wi,alpha.*phi);
shading interp;
xlabel('\omega_s (Hz)');
ylabel('\omega_i (Hz)');
title('\phi \alpha');
colorbar;
set(gca,'TickDir','out'); 

%% Integration

dz = Zmax/N;
Z = linspace(0, Zmax, dz);

trap = interp1(freqs, interp1(z, P_shift, 0), Wp).*0.5*dz + interp1(freqs, interp1(z, P_shift, Zmax), Wp).*exp(1i*delta_beta.*Zmax)*0.5*dz;
%dummy_trap = alpha.*0.5*dz + alpha.*exp(1i*delta_beta.*Zmax)*0.5*dz;

for c=1:N-1
    trap = trap + interp1(freqs, interp1(z, P_shift, dz*c), Wp).*exp(1i*delta_beta.*dz*c)*dz;
    %dummy_trap = dummy_trap + alpha.*exp(1i*delta_beta.*dz*c)*dz;
end

figure

pcolor(Ws,Wi,abs(trap));
shading interp;
xlabel('\omega_s (Hz)');
ylabel('\omega_i (Hz)');
title('Integral');
colorbar

%% Timer

elapsed_time = toc;
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);
