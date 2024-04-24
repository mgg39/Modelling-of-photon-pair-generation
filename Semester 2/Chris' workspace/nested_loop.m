clear all, close all 
clc

tic;  %Start of timer


%% Variales to tune

I1 = 10; %Number of C varialbes being tested
I2 = 10; %Number of pw variables being tested
I3 = 1; %Number of A values being tested
I4 = 1; %Number of lambda variables being tested 

C = linspace(0.7, 1.2, I1);  %Range of C values to test
pulsewidth = linspace(10, 100, I2); %Range of pw values to test
A = linspace(0.1, 0.1, I3); %Range of A values being tested
lambda = linspace(750, 750, I4)*10^-9;  %Range of lambda values to test

Purity = zeros(I1, I2, I3, I4);  %Empty 4d array for purity values

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

N = 2048; % N discretization points
T = 1500; %width of T window in ps
dt = T/N; 
t = [-T/2 : dt : T/2 - dt]'; % time domain in ps (ps determined by constants)

delta = (2 * pi / T) * [-N/2 : 1 : N/2 - 1]';  %Inverse time domain
delta = fftshift(delta);

u0=zeros(3*N, 1); %Defining an array to represent all pulses together
                  %F pulse represented by first N points, S represented by
                  %N+1 to 2N point, P represented by 2N+1 to 3N points

            %Other pulses remain at 0 for initial conditions

opts = odeset('RelTol',1e-4,'AbsTol',1e-6);  %Setting the tolerance for ode45

D_T = 14e-3;  %Damage threshold of LiNbO3, units of kJ/cm^2
w = 664e-7;    %Width of waveguide in cm
h = 330e-7;    %Height of waveguide in cm
theta = 70;    %Slant of waveguide wall in degrees

face = w*h + h^2/tan(theta); %Area of waveguide face in cm^2

load("photon_disp.mat"); %Data for singal and idler photons
lscan_photon = lamscan; neff_photon=neff;   
lscan_photon=lscan_photon*10^-6;   %Converting from um to m

load("pump_disp.mat"); lscan_pump = lamscan; neff_pump=neff; %Data for pump pulse 
lscan_pump=lscan_pump*10^-6;   %Converting from um to m         

wscan_photon=2*pi*c0./(lscan_photon);
wscan_pump=2*pi*c0./(lscan_pump); 

load('p_con_curve.mat');  %Data of p_con curve for calculations to locate where JSI forms

w_span = 2*N*10^12/(T*sqrt(2)); %Width of JSI axes so that all values of wp=ws+wi can be calculated

for J1=1:I1   %Nested for loop over values of C
    for J2=1:I2   %Nested for loop over values of pw
        for J3=1:I3   %Nested for loop over values of A
            for J4=1:I4  %Nested for loop over values of lambda
                
                fprintf("Current parameters: C = %.2f 1/cm, pw = %.1f ps, A = %.2f kW, lambda = %.2d nm \n", C(J1), pulsewidth(J2), A(J3), lambda(J4))

                w0=2*pi*c0/lambda(J4);

                delta_0 = 2*pi*c0/(2*lambda(J4)) - 2*pi*c0/(2*750*10^-9);

                ratio = 2*asech(1/2)/pulsewidth(J2);     %Finding the ratio between the desired pulsewidth and FWHM of a sech curve to scale t by
                %u0(1:N) = sech(t*ratio).*exp(-1i*delta_0*t).*sqrt(A(J3)); %*(1+1i)/sqrt(2);
                u0(1:N) = sech(t*ratio).*sqrt(A(J3)); %*(1+1i)/sqrt(2);

                E = sqrt(A(J3))*pi/(ratio*10^12);
                %Integral of a*sech(b*x) from -ininifty to infinity= a*pi/b

                DTBreach = false;
                if ((E/face) > D_T) %Checking if pulse breaches damage treshold
                    DTBreach = true;
                    fprintf("Damage threshold breached, simulation skipped. \n")
                end
                if DTBreach==false  %Only continue simulation of current parameters if damage threshold isnt breached

%% Fourier Frequency domain

                    zend = 2.885/C(J1);  %Scaling the length of the waveguide wrt C (linear coupling co-eff)
                    z = [0:zend/(N-1):zend]'; % Spacial domain in cm (cm determined by contsants)

                    [z, uhat] = ode45(@(z, uhat) CoupledPDEs(z,uhat,N,delta,Beta_f1,Beta_f2,Beta_s1,Beta_s2,Beta_p1,Beta_p2,kappa,gamma,C(J1)), z, u0, opts);
            %ode45 uses Runge-Kutta 4th and 5th order

                    u1 = uhat(:,1:N);       %Breaking apart final matrix into 3 respective pulses
                    u2 = uhat(:,N+1:2*N);
                    u3 = uhat(:,2*N+1:3*N);

%% Converting P(z, t) to P(z, w)

                    P = fft(fftshift(u3,2), [], 2);    %Defining a matrix P as the FT of u3 (after shifting u3 to t=0)
                    P_shift = fftshift(P,2);           %Shifting P so that omega=0 is centred

                    freqs = fftshift(delta);
                    freqs = freqs*1e12 + w0;  %Scaling frequency range to be centred around omega_0

%% 
                    delta_omega = 2*pi*c0*((1/(lambda(J4))) - (1/(750*10^-9)));
                    m_prime = -0.5185;

                    for I=1:5
                        x_prime = 1.3553*10^15 + delta_omega/(2*sqrt(2*(1+m_prime^2))*sin(pi/4 - abs(atan(m_prime))));
                        m_prime = spline(Xm, M, x_prime);
                    end

                    w_centre_x = 1.3553*10^15 + delta_omega/(sqrt(2*(1+m_prime^2))*sin(pi/4 - abs(atan(m_prime))));
                    w_centre_y = 1.1562*10^15 + m_prime*delta_omega/(sqrt(2*(1+m_prime^2))*sin(pi/4 - abs(atan(m_prime))));

                    wi = linspace(w_centre_x-w_span, w_centre_x+w_span, N); %Setting range of omega for idler photons
                    ws = linspace(w_centre_y-w_span, w_centre_y+w_span, N); %Setting range of omega for signal photons

                    ns=spline(wscan_photon,neff_photon,ws);  %using spline to get range of neff that follow same relation between neff_photon and wscan_photon
                    ni=spline(wscan_photon,neff_photon,wi);  

                    [Wi,Ws] = meshgrid(wi,ws);     %Converting arrays to meshgrids
                    [Ni, Ns] = meshgrid(ni, ns);   %Ws vertical, Wi horizontal
                    Wp = Ws + Wi;                  %freq_p = freq_s + freq_i

                    Np = spline(wscan_pump,neff_pump,Wp);   %using spline to get range of neff that follow same relation between neff_pump and wscan_pump

                    Beta_s = Ws.*Ns./c0;   %beta = n*w/c
                    Beta_i = Wi.*Ni./c0; 
                    Beta_p = Wp.*Np./c0;   %calculating beta values for all 3 pulses

                    delta_beta = Beta_p - Beta_s - Beta_i;  %delta_beta = beta_p - beta_s - beta_i

%% Integration

                    zend = zend/100; %converting from cm to m
                    dz = zend/N;

                    trap = interp1(freqs, interp1(z, P_shift, 0), Wp).*0.5*dz + interp1(freqs, interp1(z, P_shift, zend), Wp).*exp(1i*delta_beta.*zend)*0.5*dz;
                    %using trapezium method for integration

                    for c=1:N-1
                        trap = trap + interp1(freqs, interp1(z, P_shift, dz*c), Wp).*exp(1i*delta_beta.*dz*c)*dz; 
                    end

%% Purity
                    j = 1;

                    svdamp = svds(trap, j);
                    prob = (svdamp).^2 / ((svdamp)' * (svdamp));
                    p = sum(prob.^2);   %calculating purity of JSI by splitting it into j functions

                    diff = 1;
                    while diff>0.001  %Repeating with increasing j until difference from last iteration below threshold
                        j = j + 1;
                        svdamp = svds(trap, j);
                        prob = (svdamp).^2 / ((svdamp)' *(svdamp));
                        p = [p, sum(prob.^2);];

                        diff = p(j-1) - p(j);
                    end

                    Purity(J1, J2, J3, J4) = p(j);  %Purity of this set of parameters is the last calculated
                    fprintf("p = %.4f \n",Purity(J1, J2, J3, J4))
                end
            end
        end
    end
end

[maxPur, Index] = max(Purity(:));  %Finding the highest purity
[CIndex, PWIndex, AIndex, LIndex] = ind2sub(size(Purity), Index); %Finding the parameters corresponding to highest purity

fprintf("The highest purity is p = %.4f with the constants C = %.2f 1/cm, pulsewidth = %.1f ps, A = %.2f kW and lambda = %.2d nm", maxPur, C(CIndex), pulsewidth(PWIndex), A(AIndex), lambda(LIndex))

%% 
 
save('Purity.mat', 'Purity', 'C', 'pulsewidth', 'A', 'lambda')
fprintf("\n !!! \n MAKE SURE YOU AREN'T OVERWRITING OLD PURITY DATA THAT YOU STILL WANT \n THIS IS A WARNING BECAUSE I FORGOT TO DO SO \n !!! \n")

%% Timer

elapsed_time = toc;
mins = floor(elapsed_time/60);
secs = rem(elapsed_time, 60);
disp(['Elapsed time: ', num2str(mins), ' minutes and ', num2str(secs), ' seconds']);