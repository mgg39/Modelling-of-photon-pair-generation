clear all, close all 
clc

T = 20; % time domain width
N = 1024; % N discretization points
dt = T / N; 
t = [-T/4 : dt : 3*T/4 - dt]'; % time domain in ps (ps determined by constants)

delta = (2 * pi / T) * [-N / 2:N / 2 - 1]';  %Inverse space domain
delta = fftshift(delta);

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

max_PW = D_T*face*10^12/A;           %Maximum allowed pulsewidth for selected pulse strength without damaging crystal
max_PA = D_T*face*10^12/pulsewidth;  %Maximum allowed pulse strength for selected pulsewidth without damaging crystal

if ((A*pulsewidth*10^-12/face) > D_T)
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
figure; %%3 plots

subplot(1,3,1)         %Plotting F pulse
pcolor(t,z,abs(u1).^2)
shading interp
hold on
%plot(Beta_f1*z, z)
hold off
xlabel('t(ps)')
ylabel('x(cm)')
colorbar
ylabel(colorbar, "Pulse energy (kW)","fontsize",10,"rotation",270)
title("F")
set(gca,'TickDir','out'); 

subplot(1,3,2)         %Plotting S pulse
pcolor(t,z,abs(u2).^2)
shading interp
hold on
%plot(Beta_s1*z, z, color='w')
hold off
xlabel('t(ps)')
ylabel('x(cm)')
colorbar
ylabel(colorbar, "Pulse energy (kW)","fontsize",10,"rotation",270)
title("S")
set(gca,'TickDir','out'); 

subplot(1,3,3)          %Plotting P pulse
pcolor(t,z,abs(u3).^2)
shading interp
hold on
%plot(Beta_s1*z, z, color='w')
hold off
xlabel('t(ps)')
ylabel('x(cm)')
colorbar
ylabel(colorbar, "Pulse energy (kW)","fontsize",10,"rotation",270)
title("P")
set(gca,'TickDir','out'); 

%% Plot Fundamental Frequency vs. z
%%% Frequency Analysis

fundamental_freq_F = zeros(length(z), 1);
fundamental_freq_S = zeros(length(z), 1);
fundamental_freq_P = zeros(length(z), 1);

for i = 1:length(z)
    % Fourier transform of each pulse at the current position
    u1_hat = fft(u1(i, :));
    u2_hat = fft(u2(i, :));
    u3_hat = fft(u3(i, :));

    % Index of the maximum amplitude in the Fourier transform
    [~, index_F] = max(abs(u1_hat));
    [~, index_S] = max(abs(u2_hat));
    [~, index_P] = max(abs(u3_hat));

    % Frequency
    fundamental_freq_F(i) = abs(delta(index_F));
    fundamental_freq_S(i) = abs(delta(index_S));
    fundamental_freq_P(i) = abs(delta(index_P));
end

figure;

subplot(3, 1, 1);
plot(z, fundamental_freq_F);
xlabel('z (cm)');
ylabel('\delta F (1/ps)');
title('Fundamental Frequency \delta vs. z - F Pulse');

subplot(3, 1, 2);
plot(z, fundamental_freq_S);
xlabel('z (cm)');
ylabel('\delta S (1/ps)');
title('Fundamental Frequency vs. z - S Pulse');

subplot(3, 1, 3);
plot(z, fundamental_freq_P);
xlabel('z (cm)');
ylabel('\delta P (1/ps)');
title('Fundamental Frequency vs. z - P Pulse');