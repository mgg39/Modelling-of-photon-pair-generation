clear all, close all
clc
%% Constants

N = 1024;

Zmax = 3;      %Length of waveguide with units of cm
pw = 10;       %Width of pulse with units of ps
A = 1;         %Amplitude of pulse with units of kW

load("photon_disp.mat");
lscan_photon = lamscan; neff_photon=neff;
load("pump_disp.mat");
lscan_pump = lamscan; neff_pump=neff;

%% Initial conditions

dt = 6*pw/N;
t = [-2*pw : dt: 4*pw];

u0 = zeros(N,1);
ratio = 2*asech(1/2)/pw;     %Finding the ratio between the desired pulsewidth and FWHM of a sech curve to scale t by

for c=1:N
    u0(c)=sqrt(A)*sech(t(c)*ratio)*(1+1i)/sqrt(2);  %Sech represents laser pulses very well
                                    %Actual hiehgt of sech curve is A^2 (due to abs(u3).^2) so
                                    %is scaled by sqrt(A)
end 

omega = fft(u0);

dz = 2*Zmax/N;
z = [-Zmax : dz : Zmax-dz];
A = 1/Zmax;
H = zeros(N,1);

for c=1:N
    if z(c)>-Zmax/2 && z(c)<Zmax/2
        H(c)=A;
    end
end

p = fft(z);

%% Dispersion coefficients

beta_pump = (neff_pump./lscan_pump.*1)*2*pi;
beta_photon = (neff_photon./lscan_photon*1)*2*pi;


B = beta_pump * beta_photon';

%% Plotting

pcolor(beta_photon, beta_pump, B)
