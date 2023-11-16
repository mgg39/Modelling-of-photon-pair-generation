% https://www.youtube.com/watch?v=BkA7ncY0b7I

clear all, close all, clc

L = 100; % length domain
N = 1024; % n discretization points
dt = L / N; % Changing dx to dt
t = -L / 2:dt:L / 2 - dt; % time domain

t=t';   % AG: make it a column vector,
        % instead of row vector (ODE45 requires column vectors)

delta = (2 * pi / L) * [-N / 2:N / 2 - 1]';
delta = fftshift(delta);

% r2
Beta_f2 = 0.83e-6;
Beta_f1 = -563.3e-3;
%a = (Beta_2 / 2);    % AG: there is no need in these parameters
%b = (Beta_1);        %  it is easier and clearer to use Beta_1 and Beta_2
                      % I updated the rhspde_AG accordingly
Beta_s1 = 533.3e-3;

kappa = 6.9e5;
gamma = 100;
C = 200;

% discrete wavenumbers
delta = (2 * pi / L) * [-N / 2:N / 2 - 1]';  % AG: let's give it a proper physical name!
delta = fftshift(delta); % re-order          % these are frequency detunings delta
                                            % from the reference frequency
                                            % (central frequency of the
                                            % pump pulse)

%% Initial conditions

u0=zeros(3*N, 1); %Defining a matrix to represent all pulses together%%
for c=1:N
    u0(c, 1)=sech(t(c));  %Sech represents laser pulses very well
end
%Other pulse/s remain at 0 for initial conditions

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);

%% Fourier Frequency domain
x = [0:20/(N-1):20]'; % Replacing x
[x, uhat] = ode45(@(x, uhat) rhspde_CW(x, t, uhat, Beta_f1,Beta_f2, kappa, Beta_s1, gamma, C, N, delta), x, u0, opts);

u1 = uhat(:,1:N);
u2 = uhat(:,N+1:2*N);
u3 = uhat(:,2*N+1:3*N);

u1=ifft(u1,N,2); %Returning to spacial and temporal domains
u2=ifft(u2,N,2);
u3=ifft(u3,N,2);

%% Plot
figure;

%{
%%Waterfall
% Real part
subplot(2, 1, 1);
h_real = waterfall(t, x(1:10:end), real(u(1:10:end, :))); % Switch t and x
set(h_real, 'LineWidth', 5, 'FaceAlpha', 0.5);
colormap(flipud(jet) / 1.5);
set(gca, 'FontSize', 32);
xlabel('Time'); % Switch labels
ylabel('Space');
zlabel('Real Part');
title('Real Part of u');
set(gcf, 'Position', [1500 500 1750 1200]);

% Imaginary part
subplot(2, 1, 2);
h_imag = waterfall(t, x(1:10:end), imag(u(1:10:end, :))); % Switch t and x
set(h_imag, 'LineWidth', 5, 'FaceAlpha', 0.5);
colormap(flipud(jet) / 1.5);
set(gca, 'FontSize', 32);
xlabel('Time'); % Switch labels
ylabel('Space');
zlabel('Imaginary Part');
title('Imaginary Part of u');
%}

%{ 
%%2D standard
% Plot the real part
subplot(2, 1, 1);
plot(t, real(u), 'b');
xlabel('Time (t)');
ylabel('Real(u)');
title('Real part of u(t)');

% Plot the imaginary part
subplot(2, 1, 2);
plot(t, imag(u), 'r');
xlabel('Time (t)');
ylabel('Imag(u)');
title('Imaginary part of u(t)');

% Adjust the plot appearance
sgtitle('Real and Imaginary Parts of u(t)');

% Plot both real and imaginary parts together
plot(t, real(u), 'b', t, imag(u), 'r');
xlabel('Time (t)');
ylabel('Real and Imaginary Parts');
legend('Real(u)', 'Imag(u)');
title('Real and Imaginary Parts of u(t)');
%}

% Plot the modulus
%plot(t, abs(u).^2, 'b');
%xlabel('Time (t)');
%ylabel('(|u|^2)');
%title('|u|^2(t)');

pcolor(t,x,abs(u1).^2)
shading interp
xlabel('t(ns)')
ylabel('x(m)')
