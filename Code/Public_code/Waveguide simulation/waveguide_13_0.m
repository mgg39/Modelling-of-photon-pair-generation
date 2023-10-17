clear
% Simulation parameters
a = 0.1615;       % Width of the waveguide (m)
b = 0.08255;       % Height of the waveguide (m)
c = 3e8;        % Speed of light (m/s)
f = 1.56e9;     % Frequency (Hz)
fs = 10*f;       %Sampling Frequency(Hz)
m = 1;          % Mode index for x-direction
n = 0;          % Mode index for y-direction
length_z=10;
ending_time = 100; % Total simulation time (s)
mu=4*pi * 10^-7; % Free space permeability
eps=  8.854 * 10^-12; %Free space permittivity
B0=1;
%Sampling frequency check (Nyquest theorem)
if fs/f>=2
    disp("fs/f="+fs/f)
else
    disp("less than Nyquest limit")
    return
end
%Period Calculation
dt=1/fs;
%Calculation of wave
fc = (1/(2*pi*sqrt(mu*eps)))*sqrt((m*pi/a)^2+(n*pi/b)^2);   % Cut-off frequency (Hz)
w=2*pi*f;
% Check if the frequency is above the cut-off frequency
if f < fc
    disp(f+"<"+fc)
    return;
end

% Calculation of propagation constants
k=w*sqrt(mu*eps);
kc = sqrt((m*pi/a)^2 + (n*pi/b)^2);   % Cut-off wavenumber (1/m)
beta = sqrt(k^2-kc^2);     % Propagation constant (1/m)

% % Create meshgrid for the waveguide
 x = linspace(0, a, 100);
 y = linspace(0, b, 100);
 z = linspace(0,length_z,100);
 [X, Y, Z] = ndgrid(x, y, z);

%Define the time-stepping loop
num_steps = round(ending_time / dt);



% Calculate the electric field distribution
Ex = zeros(size(X));   % Transverse electric field component in the x-direction
Ey = zeros(size(X));   % Transverse electric field component in the y-direction
Ez = zeros(size(X));   % Transverse electric field component in the z-direction (set to zero)

%Outer forloop looping over time and plot
for t = 1:num_steps
    % Calculate the magnetic field distribution (assuming Bz component is present)
    figure(10)
    clf reset
    Ex = ((1j.*w.*mu.*n.*pi)/(kc^2*b)) * cos((m*pi/a).*X).*sin((n*pi/b).*Y).*exp(-1j*beta.*z).*exp(-1j.*w.*t.*dt) ;
    Ey = ((-1j.*w.*mu.*m.*pi)/(kc^2*a)) * sin((m*pi/a).*X).*cos((n*pi/b).*Y).*exp(-1j*beta.*z).*exp(-1j.*w.*t.*dt) ;
    
    % Plot the electric field distribution
    for i = 1:100
    surf(X(:,:,i),Y(:,:,i),Z(:,:,i),real(Ex(:,:,i))) ;
    hold on
    surf(X(:,:,i),Y(:,:,i),Z(:,:,i),real(Ey(:,:,i))) ;
    end
    hold off
    % Define custom colormap with shifted and rescaled zero value
    view(45,45)
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title(['Electric Field Distribution (TE', num2str(m), num2str(n), ' Mode)']);
    disp(t)
end