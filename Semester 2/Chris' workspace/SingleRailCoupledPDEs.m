%Equation
function duhatdt = SingleRailCoupledPDEs(z,uhat,N,delta,bf1,bf2,bs1,bs2,kappa,gamma,C)
%Inputting u rather than fft(u)
u1 = uhat(1:N); %Splitting u into seperate pulses
u2 = uhat(N+1:2*N);

fu1 = fft(u1);
fu2 = fft(u2);

duhatdt1 = ifft( (1i*bf2/2*delta.^2.*fu1.*1) + (-1i*bf1*delta.*fu1.*1) + (1i*gamma*fft(conj(u1).*u2.*exp(1i*z.*kappa).*1)) ); 
duhatdt2 = ifft( (1i*bs2/2*delta.^2.*fu2.*1) + (-1i*bs1*delta.*fu2.*1) + (1i*gamma*fft(u1.^2)*exp(-1i*z.*kappa)./2) );

duhatdt = [duhatdt1' duhatdt2']';  %Recombining pulses into one array
end