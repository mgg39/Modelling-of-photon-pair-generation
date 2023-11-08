%Equation
function duhatdt = rhspde_AG(t,uhat,z,bf1,bf2, kappa, bs1, bp1, gamma, C)
duhatdt = size(uhat);
duhatdt(1) = 1i*bf2/2*(z.^2).*uhat(1) -bf1*z.*uhat(1) +1i*gamma*exp(1i*kappa*z.)*fft((conj(duhatdt(1)).*duhatdt(2).)); %linear and diagonal
duhatdt(2) = bs1*z.*uhat(2) +1i*C*duhatdt(3) + 1i*gamma/2 *exp(-1i*kappa*z.)*fft(duhatdt(1).^2);
duhatdt(3) = -bp1*z.*uhat(3) + 1i*C*duhatdt(2);
end