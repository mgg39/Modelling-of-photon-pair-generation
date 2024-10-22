clear all, close all 
clc

tic;  %Start of timer

T = 500; % time domain width
N = 2048; % N discretization points
dt = T/N; 
t = [-T/2 : dt : T/2 - dt]'; % time domain in ps (ps determined by constants)

delta = (2 * pi / T) * [-N/2 : 1 : N/2 - 1]';  %Inverse time domain
delta = fftshift(delta);

% Constants
Beta_f2 = 0.83e-2; %Units in ps^2/cm
Beta_s2 = 0.22e-2; %Units in ps^2/cm
Beta_p2 = 2.53e-2; %Units in ps^2/cm
Beta_f1 = 563.3e-2; %Units in ps/cm
Beta_s1 = 533.3e-2; %Units in ps/cm
Beta_p1 = 0; %563.3e-2; %Test of all pulse have same beta_1

kappa =  6.9e3; %Units in 1/cm
gamma = 10^1.5; %Units in 1/(cm*sqrt(kW))
C = 1; %Units in  1/cm, C=2 corresponds to a rail seperation of x=200nm  

%% Initial conditions

u0=zeros(3*N, 1); %Defining an array to represent all pulses together
                  %F pulse represented by first N points, S represented by
                  %N+1 to 2N point, P represented by 2N+1 to 3N points

pulsewidth = 25;                    %Pulse width of laser (timeframe already scale to ps with constants)
A = 1;                              %Amplitude of laser pulse in kiloWatts (kW scaled by constants)
ratio = 2*asech(1/2)/pulsewidth;     %Finding the ratio between the desired pulsewidth and FWHM of a sech curve to scale t by

u0(1:N) = sqrt(A)*sech(t*ratio); %*(1+1i)/sqrt(2);

%Other pulses remain at 0 for initial conditions

opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
 
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
zend = 2.885/C;  %Scaling the length of the waveguide wrt C
z = [0:zend/(N-1):zend]'; % Spacial domain in cm (cm determined by contsants)

t_span = zend*Beta_f1*2; %setting a range of t to display the plots over

[z, uhat] = ode45(@(z, uhat) CoupledPDEs(z,uhat,N,delta,Beta_f1,Beta_f2,Beta_s1,Beta_s2,Beta_p1,Beta_p2,kappa,gamma,C), z, u0, opts);

u1 = uhat(:,1:N);       %Breaking apart final matrix into 3 respective pulses
u2 = uhat(:,N+1:2*N);
u3 = uhat(:,2*N+1:3*N);


%% Plot
figure;

subplot(1,3,1)         %Plotting F pulse propagation
pcolor(t,z,abs(u1).^2)
shading interp
hold on
%plot(Beta_f1*z, z, 'k-', Beta_s1*z, z, 'k--')
hold off
xlabel('t(ps)')
xlim([-t_span*3/4 t_span*5/4])
ylabel('z(cm)')
colorbar
%caxis([0 1])
ylabel(colorbar, "Pulse intensity (kW)","fontsize",10,"rotation",270)
title("F")
set(gca,'TickDir','out'); 

subplot(1,3,2)         %Plotting S pulse propagation
pcolor(t,z,abs(u2).^2)
shading interp
hold on
%plot(Beta_f1*z, z, 'w-', Beta_s1*z, z, 'w--')
hold off
xlabel('t(ps)')
xlim([-t_span*3/4 t_span*5/4])
ylabel('z(cm)')
colorbar
ylabel(colorbar, "Pulse intensity (kW)","fontsize",10,"rotation",270)
title("S")
set(gca,'TickDir','out'); 

subplot(1,3,3)          %Plotting P pulse propagation
pcolor(t,z,abs(u3).^2)
shading interp
hold on
%plot(Beta_f1*z, z, 'w-', Beta_s1*z, z, 'w--')
hold off
xlabel('t(ps)')
xlim([-t_span*3/4 t_span*5/4])
ylabel('z(cm)')
colorbar
ylabel(colorbar, "Pulse intensity (kW)","fontsize",10,"rotation",270)
title("P")
set(gca,'TickDir','out'); 

%% Converting P(z, t) to P(z, w)

%P = fft(u3, [], 2);               %Defining a matrix P as the FT of u3 (without shifting u3 to t=0)
P = fft(fftshift(u3,2), [], 2);    %Defining a matrix P as the FT of u3 (after shifting u3 to t=0)
P_shift = fftshift(P,2);           %Shifting P so that omega=0 is centred

lambda = 750*10^-9;  %Wavelength of P photons
c = 299792458;       %Speed of light
w_0 = 2*pi*c/lambda; 

freqs = fftshift(delta);
for j=1:N
    %P_shift(j,:) = P_shift(j,:)*exp(-1i*freqs(j)*T/2);  %Test for FT
end
freqs = freqs*1e12 + w_0;  %Scaling frequency range to be centred around omega_0


figure

subplot(1,3,1)
pcolor(freqs,z,abs(P_shift).^2)   %Plotting the FT(P) pulse
shading interp
hold on
xline(w_0, 'w--')   %Plotting a vertical line at w_0 to observe the offset of the frequencies
hold off
xlabel('\omega (Hz)')
xlim([2.5105e15 2.5125e15])
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Pulse intensity (kW)","fontsize",10,"rotation",270)
title("P(z, \omega)")
set(gca,'TickDir','out'); 

subplot(1,3,2)
pcolor(freqs,z,real(P_shift))   %Plotting the RE[FT(P)] pulse
shading interp
hold on
xline(w_0, 'w--')   %Plotting a vertical line at 0 to observe the offset of teh frequencies
hold off
xlabel('\omega (Hz)')
xlim([2.5105e15 2.5125e15])
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Pulse intensity (kW)","fontsize",10,"rotation",270)
title("Re[P(z, \omega)]")
set(gca,'TickDir','out'); 

subplot(1,3,3)
pcolor(freqs,z,imag(P_shift))   %Plotting the Im[FT(P)] pulse
shading interp
hold on
xline(w_0, 'w--')   %Plotting a vertical line at 0 to observe the offset of teh frequencies
hold off
xlabel('\omega (Hz)')
xlim([2.5105e15 2.5125e15])
ylabel('z (cm)')
colorbar
ylabel(colorbar, "Pulse intensity (kW)","fontsize",10,"rotation",270)
title("Im[P(z, \omega)]")
set(gca,'TickDir','out'); 


%% Peak finder

[Pmax,Idx] = max(abs(u3(:)).^2);
[PmaxRow,PmaxCol] = ind2sub(size(abs(u3).^2), Idx);   %Finding the maximum of P(z,t)
Zmax = z(PmaxRow);
Tmax = t(PmaxCol);

fprintf("For an input laser of power %.2f kW and pulsewidth %.1d ps, the Pump pulse has a maximum amplitude of %.2d kW at z = %.2d cm and t = %.1d ps\n", A, pulsewidth, Pmax, Zmax, Tmax)

[PmaxW,IdxW] = max(abs(P_shift(:)).^2);
[PmaxWRow,PmaxWCol] = ind2sub(size(abs(P_shift).^2),IdxW);  %Finding the maximum of P(z,omega)
P_z = abs(P_shift(PmaxWRow,:)).^2;    %Finding P(z=Pmax, omega)
%w_half = 2*abs(spline(P_z(1:N/2), freqs(1:N/2), PmaxW/2) - w_0);
                                                  %%Finding FWHM of P(z,
                                                  %%omega) at peak

figure
subplot(1,2,1)
plot(freqs, real(P_shift(PmaxWRow,:)))   %Plotting real and imaginary parts of P(z,omega) to observe fine structure
xlim([2.5105e15 2.5125e15])
subplot(1,2,2)
plot(freqs, imag(P_shift(PmaxWRow,:)))
xlim([2.5105e15 2.5125e15])

%% File writer

save('PumpPulse.mat', 'P_shift', 'freqs', 'z', 'P_z');   %Saving the FT pulse for import into the conservation of energy/momentum code

save('sanity_check.mat', 'u3', 'P');

%% Timer

elapsed_time = toc;
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);