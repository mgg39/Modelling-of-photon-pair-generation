% Initialize
t = linspace(0, 100, 101); % time values
h = 1; % step size

%arrays for solutions
F = zeros(1, 100);
S = zeros(1, 100);

F(1) = 1;  % Initial condition for F
gamma = 0.1;

% Loop over time steps
for c = 1:100
    % Calculate intermediate values for F using the Runge-Kutta method
    Fk1 = -gamma * S(c) * F(c);
    Fk2 = -gamma * (S(c) + h * Fk1/2) * (F(c) + h * Fk1/2);
    Fk3 = -gamma * (S(c) + h * Fk2/2) * (F(c) + h * Fk2/2);
    Fk4 = -gamma * (S(c) + h * Fk3) * (F(c) + h * Fk3);
    
    % Update F based on the weighted sum of intermediate values
    F(c+1) = F(c) + h * (Fk1 + 2 * Fk2 + 2 * Fk3 + Fk4) / 6;
    
    % Calculate intermediate values for S using the Runge-Kutta method
    Sk1 = gamma * (F(c))^2;
    Sk2 = gamma * (F(c) + h * Sk1/2)^2;
    Sk3 = gamma * (F(c) + h * Sk2/2)^2;
    Sk4 = gamma * (F(c) + h * Sk3)^2;
    
    % Update S based on the weighted sum of intermediate values
    S(c+1) = S(c) + h * (Sk1 + 2 * Sk2 + 2 * Sk3 + Sk4) / 6;
end

% Plot the results
plot(t, F, 'b-', t, S, 'r-')
