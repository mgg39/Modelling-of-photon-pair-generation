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

%% Peak finder

[Pmax,Idx] = max(abs(u3(:)).^2);
[PmaxRow,PmaxCol] = ind2sub(size(abs(u3).^2), Idx);
Zmax = z(PmaxRow);
Tmax = t(PmaxCol);

fprintf("For an input laser of power %.2f kW and pulsewidth %.1d ps, the Pump pulse has a maximum amplitude of %.2d kW at z = %.2d cm and t = %.1d ps\n", A, pulsewidth, Pmax, Zmax, Tmax)
