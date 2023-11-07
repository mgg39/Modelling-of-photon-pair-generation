clear all  %clearing variables from previous scripts
clc

gamma=1;
kappa=0.1; %Defining variables for coupled equations
C = 2;

xini=[1,0,0]; %Defining the 3 pulses to start at relative strengths of 1, 0 and 0 respectively

options = odeset('RelTol',1e-8,'AbsTol',1e-10); %No idea what this does lol

zmax=50; %Defining max z

[z,x]=ode45(@(z,x) rhs(z,x,gamma,kappa,C), [0 zmax], xini,options); %Runge_kutta method

[N, rows] = size(x); %Determining how many points were calculated duing the RK method
dt = 0.1; %defining time steps
omega0 = 2*pi/(N*dt); %calculating omega0
c_n = fft(x(:,1)); %calculating all c_n co-efficients
omega = linspace(-omega0/2, omega0/2, N); %creating a linspace in omega

figure
subplot(2,1,1)
plot(z,abs(x).^2) %Plotting 3 pulses across z
xlim([0 50])
ylim([0 1])
legend('F','S','P')
subplot(2,1,2)
plot(omega, c_n) %Plotting spectral Fourier transform of F 
xlim([-omega0./2 omega0./2])
ylim([-10 10])

%%

function  y= rhs(z,x,gamma,kappa,C) %Defining the function for the 3 bound equations
 
 y=zeros(size(x));

 y(1)=1i*gamma*conj(x(1)).*x(2).*exp(1i*kappa*z); %First pulse
 y(2)=0.5i*gamma*x(1).^2.*exp(-1i*kappa*z) - 1i*C*x(3); %Second harmonic pulse
 y(3)=-1i*C*x(2); %Pump pulse
 
end
