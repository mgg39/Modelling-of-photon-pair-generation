% This code is royally broken

function three_equations()
    % Define the time span
    tspan = [0, 10];
    
    % Define the initial conditions
    y0 = [1, 0, 0]; %[y1_initial, y2_initial, y3_initial];
    
    % Solve the system of differential equations using ode45
    [t, y] = ode45(@odeFunc, tspan, y0);
    
    % Extract the solutions for each variable
    y1 = y(:, 1);
    y2 = y(:, 2);
    y3 = y(:, 3);
    
    %{
    % Plot the results
    plot(t, y1, 'r', t, y2, 'g', t, y3, 'b');
    legend('y1', 'y2', 'y3');
    xlabel('Time');
    ylabel('Variables');
    title('Solution of Coupled Differential Equations');
    
    
    % Plot the real parts of the results
    subplot(2, 1, 1);
    plot(t, real(y1), 'r', t, real(y2), 'g', t, real(y3), 'b');
    legend('y1 (Real)', 'y2 (Real)', 'y3 (Real)');
    xlabel('Time');
    ylabel('Real Variables');
    title('Real Parts of Coupled Differential Equations');
    %}
    
    % Plot the imaginary parts of the results
    subplot(2, 1, 2);
    plot(t, imag(y1), 'r', t, imag(y2), 'g', t, imag(y3), 'b');
    legend('F', 'S', 'P'); %NOT SURE P AND F ARE THE RIGHT WAY AROUND
    xlabel('Time');
    ylabel('Imaginary Variables');
    title('Imaginary Parts of Coupled Differential Equations');
    
end


function dydt = odeFunc(t, x)
    % Function variables
    gamma = 1;
    kappa = 0.1;
    g = 1; % Coupling constant

    % Define your system of differential equations here
    dy1dt = 1i * gamma * conj(x(1)) * x(2) * exp(1i * kappa * t);
    dy2dt = 0.5i * gamma * x(1)^2 * exp(-1i * kappa * t);
    dy3dt = g * x(2);
    
    dydt = [dy1dt; dy2dt; dy3dt];
end

%main();

