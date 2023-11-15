%Equation
function duhatdt = rhspde_AG(t,uhat,z,bf1,bf2, kappa, bs1, gamma, C, N, delta)
%Inputting u rather than fft(u)
u1 = uhat(1:N); %Splitting u into seperate pulses
u2 = uhat(N+1:2*N);
u3 = uhat(2*N+1:3*N);


duhatdt1 = (1i*bf2/2*delta.^2.*u1.*1) + (1i*bf1*delta.*u1.*1) + (1i*gamma*exp(1i*z.*kappa).*fft((ifft(conj(u1)).*ifft(u2).*1)).*1); %linear and diagonal
duhatdt2 = (-1i*bs1*delta.*u2.*1) + (1i*gamma/2 *exp(1i*z.*kappa).*fft(ifft(u1).^2)) + (1i*u2.*kappa) + (1i*u3.*C);
duhatdt3 = (1i*u3.*kappa) + (1i*u2.*C);

duhatdt = [duhatdt1 duhatdt2 duhatdt3]';  %Recombining pulses into one array
end