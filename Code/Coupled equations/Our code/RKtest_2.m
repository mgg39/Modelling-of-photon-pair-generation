t = linspace(0, 100, 200);
a = 1;
F = zeros(1,200);
S = zeros(1,200);
F(1) = (1 + i)/sqrt(2);
gamma = 0.1;

for c=1:199
    Fk1 = -gamma*S(c)*conj(F(c));
    Fk2 = -gamma*(S(c) + a*Fk1/2)*conj(F(c) + a*Fk1/2);
    Fk3 = -gamma*(S(c) + a*Fk2/2)*conj(F(c) + a*Fk2/2);
    Fk4 = -gamma*(S(c) + a*Fk3)*conj(F(c) + a*Fk3);
    F(c+1) = F(c) + a*(Fk1 + 2*Fk2 + 2*Fk3 + Fk4)/6;
    Sk1 = 0.5*gamma*(F(c))^2;
    Sk2 = 0.5*gamma*(F(c) + a*Sk1/2)^2;
    Sk3 = 0.5*gamma*(F(c) + a*Sk2/2)^2;
    Sk4 = 0.5*gamma*(F(c) + a*Sk3)^2;
    S(c+1) = S(c) + a*(Sk1 + 2*Sk2 + 2*Sk3 + Sk4)/6;
end

plot(t, abs(F), 'b-', t, abs(S), 'r-')