Start_Value = 0;
End_Value = 400;
N = 800;
a = (End_Value - Start_Value)/(N-1);
z = Start_Value:a:End_Value;
F = zeros(1,N);
S = zeros(1,N);
F(1) = (1 + 1i)/sqrt(2);
gamma = 0.1;
kappa = 0.001;

for c=1:N-1
    Fk1 = 1i*gamma*S(c)*conj(F(c))*exp(1i*kappa*z(c));
    Fk2 = 1i*gamma*(S(c) + a*Fk1/2)*conj(F(c) + a*Fk1/2)*exp(1i*kappa*z(c));
    Fk3 = 1i*gamma*(S(c) + a*Fk2/2)*conj(F(c) + a*Fk2/2)*exp(1i*kappa*z(c));
    Fk4 = 1i*gamma*(S(c) + a*Fk3)*conj(F(c) + a*Fk3)*exp(1i*kappa*z(c));
    F(c+1) = F(c) + a*(Fk1 + 2*Fk2 + 2*Fk3 + Fk4)/6;
    Sk1 = 1i*0.5*gamma*(F(c))^2*exp(-1i*kappa*z(c));
    Sk2 = 1i*0.5*gamma*(F(c) + a*Sk1/2)^2*exp(-1i*kappa*z(c));
    Sk3 = 1i*0.5*gamma*(F(c) + a*Sk2/2)^2*exp(-1i*kappa*z(c));
    Sk4 = 1i*0.5*gamma*(F(c) + a*Sk3)^2*exp(-1i*kappa*z(c));
    S(c+1) = S(c) + a*(Sk1 + 2*Sk2 + 2*Sk3 + Sk4)/6;
end

plot(z, abs(F), 'b-', z, abs(S), 'r--')
xlim([Start_Value End_Value])
ylim([0 1.4])
title(['N = ',num2str(N),', a = ',num2str(a),', gamma = ',num2str(gamma), ', kappa = ',num2str(kappa)])
xlabel('z')
ylabel('Amplitude of signal')
legend({'|F|','|S|'},'Location','northwest')