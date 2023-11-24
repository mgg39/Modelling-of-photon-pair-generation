% https://www.youtube.com/watch?v=BkA7ncY0b7I

clear all, close all 
clc

L = 100; % length domain
N = 1024; % n discretization points
dt = L / N; 
t = [-L / 2:dt:L / 2 - dt]'; % time domain in ns
%t = t*10^-9; %Converting to s

delta = (2 * pi / L) * [-N / 2:N / 2 - 1]';  %Inverse space domain
delta = fftshift(delta);

% Constants
Beta_f2 = -0.83e-6; %Units in ns^2/cm
Beta_s2 = 0.22e-6; %Units in ns^2/cm
Beta_p2 = -2.53e-6; %Units in ns^2/cm
Beta_f1 = 563.3e-3; %Units in ns/cm
Beta_s1 = 533.3e-3; %Units in ns/cm
kappa = -6.9e3; %Units in 1/cm
gamma = 1*10^1.5; %Units in 1/(cm*sqrt(kW))

C = 2; %Units in  1/cm

%% Initial conditions

u0=zeros(3*N, 1); %Defining an array to represent all pulses together
                  %F pulse represented by first N points, S represented by
                  %N+1 to 2N point, P represented by 2N+1 to 3N points

pulsewidth = 10;                    %Pulse width of laser = 10ns (timeframe already scale to ns with co-efficients)
A = 0.8;                              %Amplitude of laser pulse in Watts 
FWHM = pulsewidth/(2*asech(1/2));   %Finding the full width at half maximum of the sech curve to scale u0 to the correct pulse width

for c=1:N
    u0(c)=sqrt(A)*sech(t(c)/FWHM)*(1+1i)/sqrt(2);  %Sech represents laser pulses very well
                                    %Actual hiehgt of sech curve is A^2 so
                                    %is scaled by sqrt(A)
end 
%Other pulses remain at 0 for initial conditions

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
 
%% Fourier Frequency domain
zend = 1.5;
z = [0:zend/(N-1):zend]'; % Spacial domain in cm

[z, uhat] = ode45(@(z, uhat) CoupledPDEs(z,uhat,N,delta,Beta_f1,Beta_f2,Beta_s1,Beta_s2,Beta_p2,kappa,gamma,C), z, u0, opts);

u1 = uhat(:,1:N);       %Breaking apart final matrix into 3 respective pulses
u2 = uhat(:,N+1:2*N);
u3 = uhat(:,2*N+1:3*N);


%% Plot
figure;

subplot(1,3,1)         %Plotting F pulse
pcolor(t,z,abs(u1).^2)
shading interp
xlabel('t(ns)')
ylabel('x(cm)')
colorbar
ylabel(colorbar, "Pulse energy (kW)","fontsize",10,"rotation",270)
title("F")
set(gca,'TickDir','out'); 

subplot(1,3,2)         %Plotting S pulse
pcolor(t,z,abs(u2).^2)
shading interp
xlabel('t(ns)')
ylabel('x(cm)')
colorbar
ylabel(colorbar, "Pulse energy (kW)","fontsize",10,"rotation",270)
title("S")
set(gca,'TickDir','out'); 

subplot(1,3,3)          %Plotting P pulse
pcolor(t,z,abs(u3).^2)
shading interp
xlabel('t(ns)')
ylabel('x(cm)')
colorbar
ylabel(colorbar, "Pulse energy (kW)","fontsize",10,"rotation",270)
title("P")
set(gca,'TickDir','out'); 

%% Peak finder

[P,Idx] = max(abs(u3(:)).^2);
[PmaxRow,PmaxCol] = ind2sub(size(abs(u3).^2), Idx);
Z = z(PmaxRow);
T = t(PmaxCol);

fprintf("The Pump pulse has a relative maximum of %d at z = %d cm and t = %d ns", P, Z, T)
