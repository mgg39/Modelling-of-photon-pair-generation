load('purity2.mat');

C = linspace(0.7, 1.2, 10)*100;
pw = linspace(10, 100, 10);

figure
pcolor(pw, C, Purity)
shading interp
xlabel('\tau / ps')
ylabel('C / m^{-1}')
title('A = 0.1kW, \lambda_F = 1.5\mu m')
colorbar