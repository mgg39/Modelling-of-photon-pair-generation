clear all, close all 
clc

T = 20; % time domain width
N = 1024; % N discretization points
dt = T / N; 
t = [-T/4 : dt : 3*T/4 - dt]'; % time domain in ps (ps determined by constants)

delta = (2 * pi / T) * (1:N)';

% Constants
Beta_f2 = 0.83e-2; %Units in ps^2/cm
Beta_s2 = 0.22e-2; %Units in ps^2/cm
Beta_p2 = 2.53e-2; %Units in ps^2/cm
Beta_f1 = 563.3e-2; %Units in ps/cm
Beta_s1 = 533.3e-2; %Units in ps/cm
kappa =  6.9e3; %Units in 1/cm
gamma = 1*10^1.5; %Units in 1/(cm*sqrt(kW))

C = 1; %Units in  1/cm, C=2 corresponds to a rail seperation of x=200nm  

%% Initial conditions

u0=zeros(3*N, 1); %Defining an array to represent all pulses together
                  %F pulse represented by first N points, S represented by
                  %N+1 to 2N point, P represented by 2N+1 to 3N points

pulsewidth = 1;                    %Pulse width of laser (timeframe already scale to ps with constants)
A = 1;                              %Amplitude of laser pulse in kiloWatts (kW scaled by constants)
ratio = 2*asech(1/2)/pulsewidth;     %Finding the ratio between the desired pulsewidth and FWHM of a sech curve to scale t by

for c=1:N
    u0(c)=sqrt(A)*sech(t(c)*ratio)*(1+1i)/sqrt(2);  %Sech represents laser pulses very well
                                    %Actual hiehgt of sech curve is A^2 (due to abs(u3).^2) so
                                    %is scaled by sqrt(A)
end 

%Other pulses remain at 0 for initial conditions

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
 
%% Damage threshold

D_T = 14e-3;  %Damage threshold of LiNbO3, units of kJ/cm^2
w = 664e-7;    %Width of waveguide in cm
h = 330e-7;    %Height of waveguide in cm
theta = 70;    %Slant of waveguide wall in degrees

face = w*h + h^2/tan(theta); %Area of waveguide face in cm^2

%Integral of a*sech(b*x) from -ininifty to infinity= a*pi/b
E = sqrt(A)*pi/(ratio*10^12);

max_PW = 2*asech(1/2)*D_T*face/(pi*sqrt(A)*10^-12);           %Maximum allowed pulsewidth for selected pulse strength without damaging crystal
max_PA = (D_T*face*ratio*10^12/pi)^2;  %Maximum allowed pulse strength for selected pulsewidth without damaging crystal

if ((E/face) > D_T)
    fprintf("The input laser has exceeded the damage threshold of Lithium Niobate \n")
    fprintf("Try changing the pulsewidth to be less than %.3f ps or the pulse strength to be less than %.3f kW \n", max_PW, max_PA)
    stop
end

%% Fourier Frequency domain
zend = 5*T/(8*Beta_f1);
z = [0:zend/(N-1):zend]'; % Spacial domain in cm (cm determined by contsants)

[z, uhat] = ode45(@(z, uhat) CoupledPDEs(z,uhat,N,delta,Beta_f1,Beta_f2,Beta_s1,Beta_s2,Beta_p2,kappa,gamma,C), z, u0, opts);

u1 = uhat(:,1:N);       %Breaking apart final matrix into 3 respective pulses
u2 = uhat(:,N+1:2*N);
u3 = uhat(:,2*N+1:3*N);


%% Plot
figure;

ymin = 0;
ymax = 1.1;

subplot(1,3,1)         %Plotting F pulse
pcolor(t,z,abs(u1).^2)
shading interp
hold on
%plot(Beta_f1*z, z)
hold off
xlabel('t (ps)')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Pulse power (kW)","fontsize",10,"rotation",270)
title("F")
set(gca,'TickDir','out'); 
ylim([ymin, ymax]);

subplot(1,3,2)         %Plotting S pulse
pcolor(t,z,abs(u2).^2)
shading interp
hold on
%plot(Beta_s1*z, z, color='w')
hold off
xlabel('t (ps)')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Pulse power (kW)","fontsize",10,"rotation",270)
title("S")
set(gca,'TickDir','out'); 
ylim([ymin, ymax]);

subplot(1,3,3)          %Plotting P pulse
pcolor(t,z,abs(u3).^2)
shading interp
hold on
%plot(Beta_s1*z, z, color='w')
hold off
xlabel('t (ps)')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Pulse power (kW)","fontsize",10,"rotation",270)
title("P")
set(gca,'TickDir','out'); 
ylim([ymin, ymax]);

%% Fourier Transform at each point z

U1_hat = fft(u1, N, 2); % Fourier transform of F pulse
U2_hat = fft(u2, N, 2); % Fourier transform of S pulse
U3_hat = fft(u3, N, 2); % Fourier transform of P pulse

%% Frequencies

%F
omega_F0 = 2*pi*c/1.5; %um
omegaF = omega_F0+delta/T;

% S& P
omega_SP0 = 2*omega_F0;
omega = omega_SP0+2*delta/T;

%% Wavelengths

c = physconst('LightSpeed'); % v = c

%F
wavelengthF = 2*pi*c/omegaF; %um

% S& P
wavelength = 2*pi*c/omega; %um

%% Plot Fourier Transforms

figure; % frequency

xminf = 4297;
xmaxf = 4298;
xmin = 8594;
xmax = 8596;
ymin = 0;
ymax = 1.1;

subplot(1,3,1) % Plotting Fourier transform of F pulse
pcolor(omegaF, z, abs(fftshift(U1_hat)).^2)
shading interp
xlabel('Frequency')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Field Spectrum","fontsize",10,"rotation",270)
title("F")
set(gca,'TickDir','out');
xlim([xminf, xmaxf]);
ylim([ymin, ymax]);

subplot(1,3,2) % Plotting Fourier transform of S pulse
pcolor(omega, z,abs(fftshift(U2_hat)).^2)
shading interp
xlabel('Frequency')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Field Spectrum","fontsize",10,"rotation",270)
title("S")
set(gca,'TickDir','out');
xlim([xmin, xmax]);
ylim([ymin, ymax]);

subplot(1,3,3) % Plotting Fourier transform of P pulse
pcolor(omega, z, abs(fftshift(U3_hat)).^2)
shading interp
xlabel('Frequency')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Field Spectrum","fontsize",10,"rotation",270)
title("P")
set(gca,'TickDir','out');
xlim([xmin, xmax]);
ylim([ymin, ymax]);

figure; %wavelengths

% Dynamically set x-axis limits based on the calculated wavelengths
xminf = min(wavelengthF(:));
xmaxf = max(wavelengthF(:));
xmin = min(wavelength(:));
xmax = max(wavelength(:));

subplot(1,3,1) % Plotting Fourier transform of F pulse
pcolor(wavelengthF, z, abs(fftshift(U1_hat)).^2) % wavelength = v / f
shading interp
xlabel('Wavelength')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Field Spectrum","fontsize",10,"rotation",270)
title("F")
set(gca,'TickDir','out');
xlim([xminf, xmaxf]);

subplot(1,3,2) % Plotting Fourier transform of S pulse
pcolor(wavelength, z,abs(fftshift(U2_hat)).^2)
shading interp
xlabel('Wavelength')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Field Spectrum","fontsize",10,"rotation",270)
title("S")
set(gca,'TickDir','out');
xlim([xmin, xmax]);

subplot(1,3,3) % Plotting Fourier transform of P pulse
pcolor(wavelength, z, abs(fftshift(U3_hat)).^2)
shading interp
xlabel('Wavelength')
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Field Spectrum","fontsize",10,"rotation",270)
title("P")
set(gca,'TickDir','out');
xlim([xmin, xmax]);