t = linspace(0, 100, 101);
h = 1;
F = zeros(1,100);
S = zeros(1,100);
F(1) = 1;
gamma = 0.1;

for c=1:100
    Fk1 = -gamma*S(c)*F(c);
    Fk2 = -gamma*(S(c) + h*Fk1/2)*(F(c) + h*Fk1/2);
    Fk3 = -gamma*(S(c) + h*Fk2/2)*(F(c) + h*Fk2/2);
    Fk4 = -gamma*(S(c) + h*Fk3)*(F(c) + h*Fk3);
    F(c+1) = F(c) + h*(Fk1 + 2*Fk2 + 2*Fk3 + Fk4)/6;
    Sk1 = gamma*(F(c))^2;
    Sk2 = gamma*(F(c) + h*Sk1/2)^2;
    Sk3 = gamma*(F(c) + h*Sk2/2)^2;
    Sk4 = gamma*(F(c) + h*Sk3)^2;
    S(c+1) = S(c) + h*(Sk1 + 2*Sk2 + 2*Sk3 + Sk4)/6;
end

plot(t, F, 'b-', t, S, 'r-')