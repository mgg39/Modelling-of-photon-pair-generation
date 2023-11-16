%Equation
function duhatdt = rhspde_CW(z, t, uhat, bf1,bf2, kappa, bs1, gamma, C, N, delta)
%Inputting u rather than fft(u)
u1 = uhat(1:N); %Splitting u into seperate pulses
u2 = uhat(N+1:2*N);
u3 = uhat(2*N+1:3*N);

fu1 = fft(u1);
fu2 = fft(u2);
fu3 = fft(u3);

duhatdt1 = (1i*bf2/2*delta.^2.*fu1.*1) + (1i*bf1*delta.*fu1.*1); %+ (1i*gamma*exp(1i*z.*kappa).*fft((conj(u1).*u2.*1)).*1); %linear and diagonal
duhatdt2 = (-1i*bs1*delta.*fu2.*1) + (1i*fu2.*kappa) + (1i*fu3.*C); % + (1i*gamma/2 *exp(1i*z.*kappa).*fft(u1.^2));
duhatdt3 = (1i*fu3.*kappa) + (1i*fu2.*C);

duhatdt = [duhatdt1' duhatdt2' duhatdt3']';  %Recombining pulses into one array
end