%Equation
function duhatdt = rhspde_AG(t,uhat,z,bf1,bf2, kappa, bs1, bp1, gamma, C)
%Inputting u rather than fft(u)
u1 = uhat(:,1); %Splitting u into seperate pulses
u2 = uhat(:,2);
u3 = uhat(:,3)
ftu1 = fft(u1); %Calling ftt on each pulse
ftu2 = fft(u2);
ftu3 = fft(u3);

duhatdt1 = bf2/2*ftu2.^2 -ftu1.*bf1 + 1i*gamma*exp(1i*z.*kappa)*fft((conj(u1).*u2.*1)); %linear and diagonal
duhatdt2 = ftu2.*bs1 + 1i*gamma/2 *exp(1i*z.*kappa)*fft(u1.^2) + 1i*u2.*kappa +1i*u3.*C;
duhatdt3 = 1i*u3.*kappa + ftu3.*bp1 + 1i*u2.*C;

duhatdt = [duhatdt1 duhatdt2 duhatdt3];
end