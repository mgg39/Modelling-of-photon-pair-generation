% https://www.youtube.com/watch?v=BkA7ncY0b7I

clear all, close all, clc

L = 100; % length domain
N = 1024; % n discretization points
dt = L / N; 
t = [-L / 2:dt:L / 2 - dt]'; % time domain

delta = (2 * pi / L) * [-N / 2:N / 2 - 1]';  %Inverse space domain
delta = fftshift(delta);

% Constants
Beta_f2 = 0.83e-6;
Beta_f1 = 563.3e-3;
Beta_s1 = 533.3e-3;
kappa = 6.9e5;
gamma = 100;

C = 50;

%% Initial conditions

u0=zeros(3*N, 1); %Defining a matrix to represent all pulses together
                  %F pulse represented by first N points, S represented by
                  %N+1 to 2N point, P represented by 2N+1 to 3N points

for c=1:N
    u0(c, 1)=sech(t(c));  %Sech represents laser pulses very well
end
%Other pulses remain at 0 for initial conditions

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
 
%% Fourier Frequency domain
z = [0:20/(N-1):20]'; % Spacial domain
[z, uhat] = ode45(@(z, uhat) rhspde_CW(z,uhat, Beta_f1,Beta_f2, kappa, Beta_s1, gamma, C, N, delta), z, u0, opts);

u1 = uhat(:,1:N);       %Breaking apart final matrix into 3 respective pulses
u2 = uhat(:,N+1:2*N);
u3 = uhat(:,2*N+1:3*N);


%% Plot
figure;

subplot(1,3,1)         %Plotting F pulse
pcolor(t,z,abs(u1).^2)
shading interp
xlabel('t(ns)')
ylabel('x(m)')
colorbar
title("F")

subplot(1,3,2)         %Plotting S pulse
pcolor(t,z,abs(u2).^2)
shading interp
xlabel('t(ns)')
ylabel('x(m)')
colorbar
title("S")

subplot(1,3,3)          %Plotting P pulse
pcolor(t,z,abs(u3).^2)
shading interp
xlabel('t(ns)')
ylabel('x(m)')
colorbar
title("P")

