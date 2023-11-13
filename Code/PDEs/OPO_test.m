%% This code solves the set of coupled equation describing OPO process in chi2
%   with no dispersion in pump field (omega0)
% and only linear dispersion in signal+idler fields (~omega0/2)
% [but it is easy to modify it to account for any arbitrary dispersion!
%  simply modify the definitions of betap, beta1, and beta2]
%  
%               d_z P=-i A1 A2                   (1)
%               d_z A1=-vd_t A_1-iP conj(A2)     (2)
%               d_z A2=+vd_t A_2-iP conj(A1)     (3)
% using two methods:
%   (1) standard split-step
%   (2) a modified UPPE (Universal Pulse Propagation Equation - see DOI:
%   10.1140/epjst/e2011-01503-3)
%
% For method (2) we first re-write Eqs. (1-3) in frequency domain:
%
%              d_z Pw=i betap*Pw - i FFT[A1*A2]
%              d_z A1w=i beta1*A1w - i FFT[P*conj(A2)]
%              d_z A2w=i beta2*A2w - i FFT[P*conj(A1)]
%
%              where Aw=FFT[A] and 
%              betap=0
%              beta1,2=+-v omega
%
% Next, we change to the rotating frame by introducing new variable
%               Psi=Aw*exp(-i beta z),   Aw= Psi*exp(i beta z)
% Psi satisfies:
%             d_z PsiP= -i FFT[A1 A2]*exp(-i betap z)
%             d_z Psi_A1= - i FFT[P*conj(A2)]*exp(-i beta1 z)
%             d_z Psi_A2= - i FFT[P*conj(A1)]*exp(-i beta2 z)
%
%            with A=IFFT[Aw]=IFFT[Psi*exp(i beta z)]
%


clear all
clc

%% define grid 

N=1024;

dx=0.01;
x=[-N/2:N/2-1]'*dx;

dk=2*pi/N/dx;
k=fftshift([-N/2:N/2-1]'*dk);

L=30;    % propagation distance
steps=1000;  % number of steps to take along z 

dz=L/steps;

%% define dispersion
v=0.1;

beta1=v*k;
beta2=-v*k;
betap=zeros(size(k));

%% define ini condition

P=sech(x*10);
%A1=zeros(size(x));
%A2=zeros(size(x));

A1=1e-5*sech(x*10);
A2=1e-5*sech(x*10);

%% propagate split-step

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

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 psi=vertcat(P,A1,A2);

 sol=ode45(@(z,psi) rhs(z,psi,N), [0 dz], psi);

 psi=deval(sol,dz);

 P=psi(1:N);
 A1=psi(N+1:2*N);
 A2=psi(2*N+1:3*N);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 Pout(j+1,:)=P;
 A1out(j+1,:)=A1;
 A2out(j+1,:)=A2;

end

tsplit=toc

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

%% now try uppe

% reset initial state
P=Pout(1,:);
A1=A1out(1,:);
A2=A2out(1,:);

psi=vertcat(shiftdim(fft(P),1),shiftdim(fft(A1),1),shiftdim(fft(A2),1));

tic

sol=ode45(@(z,psi) rhs_uppe(z,psi,N,betap,beta1,beta2), [0 L], psi);

psiout=shiftdim(deval(sol,z),1);

tuppe=toc

pspout=psiout(:,1:N);
a1spout=psiout(:,N+1:2*N);
a2spout=psiout(:,2*N+1:3*N);

[BP,Z]=meshgrid(betap,shiftdim(z,1));

Pout_uppe=ifft(pspout.*exp(1i*(BP.*Z)),N,2);

[B1,Z]=meshgrid(beta1,shiftdim(z,1));

A1out_uppe=ifft(a1spout.*exp(1i*(B1.*Z)),N,2);

[B2,Z]=meshgrid(beta2,shiftdim(z,1));

A2out_uppe=ifft(a2spout.*exp(1i*(B2.*Z)),N,2);

pcolor(x,z,abs(Pout)) 
shading interp

%a=ifft(psiout.*exp(1i*(B.*Z)),N,2);  % psiout is a Zn x N matrix, where Zn is the  
                  % number of points along z. We need to peform IFT in the
                  % second dimension of this matrix


%%

function  y= rhs(z,x,N)

 p=-1i*x(N+1:2*N).*x(2*N+1:3*N);
 a1=-1i*x(1:N).*conj(x(2*N+1:3*N));
 a2=-1i*x(1:N).*conj(x(N+1:2*N));

 y=vertcat(p,a1,a2);

end

