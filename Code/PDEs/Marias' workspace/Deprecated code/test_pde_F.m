% i dFdz = (Beta_2 / 2) ddFddt

clear all
clc

Beta_2 = 0.83;

xini = [1, 0, 0];

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

zmax = 50;

[z, x] = ode45(@(z, x) pde(Beta_2, x), [0 zmax], xini, options);

figure
plot(z, abs(x).^2)
xlim([0 50])
ylim([0 1])
legend('F')

function y = pde(Beta_2, x)

y = zeros(size(x));

y(1) = (1/1i) * (-Beta_2/2) * (ddy(1)/ddt);

end
