clear all; close all; 
clc;

L = 20; % length domain
N = 1024; % N discretization points
dt = L / N; 
t = [-L / 2:dt:L / 2 - dt]'; % time domain in ns (ns determined by constants)

delta = (2 * pi / L) * [-N / 2:N / 2 - 1]';  %Inverse space domain
delta = fftshift(delta);

% Constants
Beta_f2 = -0.83e-8; %Units in ns^2/cm
Beta_s2 = 0.22e-8; %Units in ns^2/cm
Beta_p2 = -2.53e-8; %Units in ns^2/cm
Beta_f1 = 563.3e-5; %Units in ns/cm
Beta_s1 = 533.3e-5; %Units in ns/cm
kappa = -6.9e3; %Units in 1/cm
gamma = 1*10^1.5; %Units in 1/(cm*sqrt(kW))

C = 2; %Units in  1/cm

%% Initial conditions

u0=zeros(3*N, 1); %Defining an array to represent all pulses together
                  %F pulse represented by first N points, S represented by
                  %N+1 to 2N point, P represented by 2N+1 to 3N points

pulsewidth = 1;  %Pulse width of laser (timeframe already scale to ns with co-efficients)
A = 0.04;         %Amplitude of laser pulse in kiloWatts 
ratio = 2*asech(1/2)/pulsewidth; %Finding the ratio between the desired pulsewidth and FWHM of a sech curve to scale t by

for c=1:N
    u0(c)=sqrt(A)*sech(t(c)*ratio)*(1+1i)/sqrt(2);  %Sech represents laser pulses very well
    %Actual height of sech curve is A^2 (due to abs(u3).^2) so is scaled by sqrt(A)
end 

%Other pulses remain at 0 for initial conditions

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
 
%% Damage threshold

D_T = 14e-3;  %Damage threshold of LiNbO3, units of kJ/cm^2
w = 664e-7;    %Width of waveguide in cm
h = 330e-7;    %Height of waveguide in cm
theta = 70;    %Slant of waveguide wall in degrees

face = w*h + h^2/tan(theta); %Area of waveguide face in cm^2

max_PW = D_T*face*10^9/A;           %Maximum allowed pulsewidth for selected pulse strength without damaging crystal
max_PA = D_T*face*10^9/pulsewidth;  %Maximum allowed pulse strength for selected pulsewidth without damaging crystal

if ((A*pulsewidth*10^-9/face) > D_T)
    fprintf("The input laser has exceeded the damage threshold of Lithium Niobate \n")
    fprintf("Try changing the pulsewidth to be less than %.3f ns or the pulse strength to be less than %.3f kW \n", max_PW, max_PA)
    stop
end

%% Fourier Frequency domain
zend = 1.5;
z = [0:zend/(N-1):zend]'; % Spacial domain in cm (cm determined by contsants)

[z, uhat] = ode45(@(z, uhat) CoupledPDEs(z,uhat,N,delta,Beta_f1,Beta_f2,Beta_s1,Beta_s2,Beta_p2,kappa,gamma,C), z, u0, opts);

u1 = uhat(:,1:N);       %Breaking apart final matrix into 3 respective pulses
u2 = uhat(:,N+1:2*N);
u3 = uhat(:,2*N+1:3*N);

%%%%% Dispersion 

% Dispersion relation (assuming linear dispersion for simplicity)
k = linspace(-pi/dt, pi/dt, N);  % Wavenumber range

omega = sqrt(Beta_f1 + Beta_f2 * k.^2) + sqrt(Beta_s1 + Beta_s2 * k.^2) + sqrt(Beta_p2 * k.^2); 
omega = abs(omega)

% Time values for the dispersion curve
time_values = 1 ./ omega;

figure;

% Plot dispersion curve
subplot(2, 1, 1);
plot(k,omega);
xlabel('Wavenumber');
ylabel('Magnitude of Angular frequency ');
title('Dispersion Curve');

% Plot dispersion curve against time
subplot(2, 1, 2);
plot(time_values, omega);
xlabel('Time (ns)');
ylabel('Magnitude of Angular frequency ');
title('Dispersion Curve against Time');

