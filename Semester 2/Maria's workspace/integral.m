clear all, close all
clc

% Define the parameters
B = 1; % Value of B
O = 2; % Value of O

% Define the function f(z) = p(exp(iBOz))
f = @(z) p(exp(1i * B * O * z));

% Define the range of z values
z_values = linspace(-10, 10, 1000); % Adjust range and number of points as needed

% Compute the integral of f(z) with respect to z
integral_values = zeros(size(z_values));
for i = 1:length(z_values)
    integral_values(i) = compute_integral(f, 0, z_values(i)); % Change integration limits as needed
end

% Plot the real and imaginary parts of the integral
figure;
plot(z_values, real(integral_values), 'b', z_values, imag(integral_values), 'r');
xlabel('z');
ylabel('Integral Value');
title('Plot of Integral of f(z) = p(exp(iBOz))');
legend('Real Part', 'Imaginary Part');
grid on;

% Function to compute integral (renamed to avoid conflict)
function integral_value = compute_integral(f, a, b)
    integral_value = integral(f, a, b);
end
