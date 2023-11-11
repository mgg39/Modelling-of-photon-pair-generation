%Equation
function duhatdt = rhspde_AG(t,uhat,z,bf1,bf2, kappa, bs1, bp1, gamma, C)
u1 = uhat(:,1);
u2 = uhat(:,2);
ftu1 = fft(u1);
ftu2 = fft(u2);
duhatdt = size(uhat);
duhatdt(:,1) = bf2/2*ftu2.^2. -bf1*z.*ftu1. + 1i*gamma*exp(1i*kappa*z.)*fft((conj(u1).*u2.)); %linear and diagonal
duhatdt(:,2) = bs1*z.*ftu2. + 1i*gamma/2 *exp(-1i*kappa*z.)*fft(u1.^2) + 1i*kappa*u2.;
end