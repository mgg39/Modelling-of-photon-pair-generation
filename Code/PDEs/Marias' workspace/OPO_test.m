clear all
clc

%% define grid - don't touch for now

N=1024;

dx=0.01;
x=[-N/2:N/2-1]'*dx;

dk=2*pi/N/dx;
k=fftshift([-N/2:N/2-1]'*dk);

L=30;    % propagation distance
steps=1000;  % number of steps to take along z 

dz=L/steps;

%% define dispersion - unsure
v=0.1;

beta1=v*k;
beta2=-v*k;
betap=zeros(size(k));

%% define ini condition - pump: sech n - gaussian

P=sech(x*10);
%A1=zeros(size(x));
%A2=zeros(size(x));

A1=1e-5*sech(x*10);
A2=1e-5*sech(x*10);

%% propagate split-step - unsure

Pout=zeros(steps+1,N);
A1out=zeros(steps+1,N);
A2out=zeros(steps+1,N);

Pout(1,:)=P;
A1out(1,:)=A1;
A2out(1,:)=A2;

prop1=exp(1i*beta1*dz);
prop2=exp(1i*beta2*dz);
propp=exp(1i*betap*dz);


tic

for j=1:steps

 P=ifft(propp.*fft(P));
 A1=ifft(prop1.*fft(A1));
 A2=ifft(prop2.*fft(A2));

 %%%%%%%%%%%%CHECK%%%%%%%%%%%%%%%%%%
 psi=vertcat(P,A1,A2);

 sol=ode45(@(z,psi) rhs(z,psi,N), [0 dz], psi);

 psi=deval(sol,dz);

 P=psi(1:N);
 A1=psi(N+1:2*N);
 A2=psi(2*N+1:3*N);

 %%%%%%%%%%%%CHECK%%%%%%%%%%%%%%%%%%
 
 Pout(j+1,:)=P;
 A1out(j+1,:)=A1;
 A2out(j+1,:)=A2;

end

tsplit=toc

%%%%%%%%%%%%%%%%%%PLOTS%%%%%%%%%%%%%%%%%%%%%%%%
%% plots

z=[0:steps]*dz;

figure              % plot the field in time domain
pcolor(x,z,abs(Pout)) 
shading interp
xlabel('t')
ylabel('z')
title('pump')

figure            % plot the spectrum
pcolor(fftshift(k),z,abs(fftshift(fft(Pout,N,2),2)))
shading interp
xlabel('k')
ylabel('z')
title('pump spectrum')

figure              % plot the field in time domain
pcolor(x,z,abs(A1out)) 
shading interp
xlabel('t')
ylabel('z')
title('A1')

figure            % plot the spectrum
pcolor(fftshift(k),z,abs(fftshift(fft(A1out,N,2),2)))
shading interp
xlabel('k')
ylabel('z')
title('A1 spectrum')

figure              % plot the field in time domain
pcolor(x,z,abs(A2out)) 
shading interp
xlabel('t')
ylabel('z')
title('A2')

figure            % plot the spectrum
pcolor(fftshift(k),z,abs(fftshift(fft(A2out,N,2),2)))
shading interp
xlabel('k')
ylabel('z')
title('A2 spectrum')

ppow=sum(abs(Pout).^2,2)*dx;
a1pow=sum(abs(A1out).^2,2)*dx;
a2pow=sum(abs(A2out).^2,2)*dx;

figure
plot(z,ppow,z,a1pow,z,a2pow)
xlabel('z')
ylabel('Energy')

%%%%%%%%%%%%%%%%%%PLOTS%%%%%%%%%%%%%%%%%%%%%%%%


function  y= rhs(z,x,N)

 p=-1i*x(N+1:2*N).*x(2*N+1:3*N);
 a1=-1i*x(1:N).*conj(x(2*N+1:3*N));
 a2=-1i*x(1:N).*conj(x(N+1:2*N));

 y=vertcat(p,a1,a2);

end