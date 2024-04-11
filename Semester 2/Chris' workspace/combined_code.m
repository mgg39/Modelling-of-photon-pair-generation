clear all, close all 
clc

tic;  %Start of timer


%% Variales to tune

pulsewidth = 25;
C = 1; %Units in  1/cm, C=2 corresponds to a rail seperation of x=200nm  
A = 1; %Amplitude of laser pulse in kiloWatts (kW scaled by constants)
lambda = 750*10^-9;  %Wavelength of P photons

PLOT = true;

%% Constants
Beta_f2 = 0.83e-2; %Units in ps^2/cm
Beta_s2 = 0.22e-2; %Units in ps^2/cm
Beta_p2 = 2.53e-2; %Units in ps^2/cm
Beta_f1 = 563.3e-2; %Units in ps/cm
Beta_s1 = 533.3e-2; %Units in ps/cm
Beta_p1 = 0; %563.3e-2; %Test of all pulse have same beta_1

kappa =  6.9e3; %Units in 1/cm
gamma = 10^1.5; %Units in 1/(cm*sqrt(kW))

c0 = 299792458; %Speed of light

w0=2*pi*c0/lambda;

T = 1500;
N = 2048; % N discretization points
dt = T/N; 
t = [-T/2 : dt : T/2 - dt]'; % time domain in ps (ps determined by constants)

delta = (2 * pi / T) * [-N/2 : 1 : N/2 - 1]';  %Inverse time domain
delta = fftshift(delta);

%% Initial conditions

u0=zeros(3*N, 1); %Defining an array to represent all pulses together
                  %F pulse represented by first N points, S represented by
                  %N+1 to 2N point, P represented by 2N+1 to 3N points

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
zend = 2.885/C;  %Scaling the length of the waveguide wrt C (linear coupling co-eff)
z = [0:zend/(N-1):zend]'; % Spacial domain in cm (cm determined by contsants)

t_span = zend*Beta_f1*4; %setting a range of t to display the plots over

[z, uhat] = ode45(@(z, uhat) CoupledPDEs(z,uhat,N,delta,Beta_f1,Beta_f2,Beta_s1,Beta_s2,Beta_p1,Beta_p2,kappa,gamma,C), z, u0, opts);

u1 = uhat(:,1:N);       %Breaking apart final matrix into 3 respective pulses
u2 = uhat(:,N+1:2*N);
u3 = uhat(:,2*N+1:3*N);


%% Plot

if (PLOT == true)
    figure;

    subplot(1,3,1)         %Plotting F pulse propagation
    pcolor(t,z,abs(u1).^2)
    shading interp
    hold on
    %plot(Beta_f1*z, z, 'k-', Beta_s1*z, z, 'k--')
    hold off
    xlabel('t(ps)', "fontSize", 12)
    xlim([-t_span*3/4 t_span*5/4])
    ylabel('z(cm)', "fontSize", 12)
    ylabel(colorbar, "Pulse intensity (kW)","fontsize",12,"rotation",270,"position",[3.5 0.5])
    title("F", "fontsize", 15)
    set(gca,'TickDir','out'); 

    subplot(1,3,2)         %Plotting S pulse propagation
    pcolor(t,z,abs(u2).^2)
    shading interp
    hold on
    %plot(Beta_f1*z, z, 'w-', Beta_s1*z, z, 'w--')
    hold off
    xlabel('t(ps)', "fontSize", 12)
    xlim([-t_span*3/4 t_span*5/4])
    ylabel('z(cm)', "fontSize", 12)
    ylabel(colorbar, "Pulse intensity (kW)","fontsize",12,"rotation",270,"position",[3.5 0.5])
    title("S", "fontsize", 15)
    set(gca,'TickDir','out'); 

    subplot(1,3,3)          %Plotting P pulse propagation
    pcolor(t,z,abs(u3).^2)
    shading interp
    hold on
    %plot(Beta_f1*z, z, 'w-', Beta_s1*z, z, 'w--')
    hold off
    xlabel('t(ps)', "fontSize", 12)
    xlim([-t_span*3/4 t_span*5/4])
    ylabel('z(cm)', "fontSize", 12)
    ylabel(colorbar, "Pulse intensity (kW)","fontsize",12,"rotation",270,"position",[3.5 0.5])
    title("P", "fontsize", 15)
    set(gca,'TickDir','out'); 
end

%% Converting P(z, t) to P(z, w)

P = fft(fftshift(u3,2), [], 2);    %Defining a matrix P as the FT of u3 (after shifting u3 to t=0)
P_shift = fftshift(P,2);           %Shifting P so that omega=0 is centred

freqs = fftshift(delta);
freqs = freqs*1e12 + w0;  %Scaling frequency range to be centred around omega_0

freqs_range = (max(freqs) - min(freqs))/2;

if (PLOT == true)
    figure

    subplot(1,3,1)
    pcolor(freqs,z,abs(P_shift).^2)   %Plotting the FT(P) pulse
    shading interp
    hold on
    %xline(w0, 'w--')   %Plotting a vertical line at w_0 to observe the offset of the frequencies
    hold off
    xlabel('\omega (Hz)', "fontsize", 12)
    xlim([w0-freqs_range w0+freqs_range])
    ylabel('z (cm)', "fontsize", 12)
    ylabel(colorbar, "Pulse intensity (kW)","fontsize",12,"rotation",270,"position",[3.5 0.5])
    title("P(z, \omega)", "fontsize", 15)
    set(gca,'TickDir','out'); 

    subplot(1,3,2)
    pcolor(freqs,z,real(P_shift))   %Plotting the RE[FT(P)] pulse
    shading interp
    hold on
    %xline(w0, 'w--')   %Plotting a vertical line at 0 to observe the offset of teh frequencies
    hold off
    xlabel('\omega (Hz)', "fontsize", 12)
    xlim([w0-freqs_range w0+freqs_range])
    ylabel('z (cm)', "fontsize", 12)
    ylabel(colorbar, "Pulse intensity (kW)","fontsize",12,"rotation",270,"position",[3.5 0.5])
    title("Re[P(z, \omega)]", "fontsize", 15)
    set(gca,'TickDir','out'); 

    subplot(1,3,3)
    pcolor(freqs,z,imag(P_shift))   %Plotting the Im[FT(P)] pulse
    shading interp
    hold on
    %xline(w0, 'w--')   %Plotting a vertical line at 0 to observe the offset of teh frequencies
    hold off
    xlabel('\omega (Hz)', "fontsize", 12)
    xlim([w0-freqs_range w0+freqs_range])
    ylabel('z (cm)', "fontsize", 12)
    ylabel(colorbar, "Pulse intensity (kW)","fontsize",12,"rotation",270,"position",[3.5 0.5])
    title("Im[P(z, \omega)]", "fontsize", 15)
    set(gca,'TickDir','out'); 
end

%% Peak finder

[Pmax,Idx] = max(abs(u3(:)).^2);
[PmaxRow,PmaxCol] = ind2sub(size(abs(u3).^2), Idx);   %Finding the maximum of P(z,t)
Zmax = z(PmaxRow);
Tmax = t(PmaxCol);

fprintf("For an input laser of power %.2f kW and pulsewidth %.1d ps, the Pump pulse has a maximum amplitude of %.2d kW at z = %.2d cm and t = %.1d ps\n", A, pulsewidth, Pmax, Zmax, Tmax)

[PmaxW,IdxW] = max(abs(P_shift(:)).^2);
[PmaxWRow,PmaxWCol] = ind2sub(size(abs(P_shift).^2),IdxW);  %Finding the maximum of P(z,omega)

if (PLOT == true)
    figure
    subplot(1,2,1)
    plot(freqs, real(P_shift(PmaxWRow,:)))   %Plotting real and imaginary parts of P(z,omega) to observe fine structure
    xlim([w0-freqs_range w0+freqs_range])
    subplot(1,2,2)
    plot(freqs, imag(P_shift(PmaxWRow,:)))
    xlim([w0-freqs_range w0+freqs_range])
end

%% Loading photon data

load("photon_disp.mat"); %Data for i and s pulses
lscan_photon = lamscan; neff_photon=neff;   
lscan_photon=lscan_photon*10^-6;   %Converting from um to m

load("pump_disp.mat"); lscan_pump = lamscan; neff_pump=neff; %Data for pump pulse 
lscan_pump=lscan_pump*10^-6;   %Converting from um to m

%% Creating meshgrids

load('p_con_curve.mat');

m_prime = -0.5185;
delta_omega = 2*pi*c0*((1/(lambda)) - (1/(750*10^-9)));

for I=1:5
    x_prime = 1.3553*10^15 + delta_omega/(2*sqrt(2*(1+m_prime^2))*sin(pi/4 - abs(atan(m_prime))));
    m_prime = spline(Xm, M, x_prime);
end

w_centre_x = 1.3553*10^15 + delta_omega/(sqrt(2*(1+m_prime^2))*sin(pi/4 - abs(atan(m_prime))));
w_centre_y = 1.1562*10^15 + m_prime*delta_omega/(sqrt(2*(1+m_prime^2))*sin(pi/4 - abs(atan(m_prime))));

w_span = 2*N*10^12/(T*sqrt(2));

wi = linspace(w_centre_x-w_span, w_centre_x+w_span, N); %Setting range of omega for idler photons
ws = linspace(w_centre_y-w_span, w_centre_y+w_span, N); %Setting range of omega for signal photons

wscan_photon=2*pi*c0./(lscan_photon);

ns=spline(wscan_photon,neff_photon,ws);  %using spline to get range of neff that follow same relation between neff_photon and wscan_photon
ni=spline(wscan_photon,neff_photon,wi);  

wscan_pump=2*pi*c0./(lscan_pump); 

[Wi,Ws] = meshgrid(wi,ws);     %Converting arrays to meshgrids
[Ni, Ns] = meshgrid(ni, ns);   %Ws vertical, Wi horizontal
Wp = Ws + Wi;                  %freq_p = freq_s + freq_i

Np = spline(wscan_pump,neff_pump,Wp);   %using spline to get range of neff that follow same relation between neff_pump and wscan_pump

%% Calculating betas

Beta_s = Ws.*Ns./c0;   %beta = n*w/c
Beta_i = Wi.*Ni./c0; 
Beta_p = Wp.*Np./c0;   %calculating beta values for all 3 pulses

delta_beta = Beta_p - Beta_s - Beta_i;  %delta_beta = beta_p - beta_s - beta_i

%% Integration

zend = zend/100;
dz = zend/N;

trap = interp1(freqs, interp1(z, P_shift, 0), Wp).*0.5*dz + interp1(freqs, interp1(z, P_shift, zend), Wp).*exp(1i*delta_beta.*zend)*0.5*dz;

for c=1:N-1
    trap = trap + interp1(freqs, interp1(z, P_shift, dz*c), Wp).*exp(1i*delta_beta.*dz*c)*dz; 
end

%% 


if (PLOT==true)
    figure
    pcolor(Wi,Ws,abs(trap).^2);
    shading interp;
    hold on
    %plot(ws, w0-ws, 'w--')
    hold off
    xlabel('\omega_s (Hz)', "fontsize", 12);
    ylabel('\omega_i (Hz)', "fontsize", 12);
    title('JSI', "fontsize", 15);
end

%% Purity

j = 1;

svdamp = svds(trap, j);
prob = (svdamp).^2 / ((svdamp)' * (svdamp));
p = sum(prob.^2);   % purity

diff = 1;
while diff>0.001
    j = j + 1;
    svdamp = svds(trap, j);
    prob = (svdamp).^2 / ((svdamp)' *(svdamp));
    p = [p, sum(prob.^2);];

    diff = p(j-1) - p(j);
end

figure
plot(p)
title(['P = ', num2str(p(j)), ', j = ', num2str(j)])
xlabel('Number of funcitons decomposed into (j)')
ylabel('Purity')
xlim([0 j])
ylim([0 1])

%% Timer

elapsed_time = toc;
mins = floor(elapsed_time/60);
secs = rem(elapsed_time, 60);
disp(['Elapsed time: ', num2str(mins), ' minutes and ', num2str(secs), ' seconds']);