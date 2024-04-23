load('Purity.mat');

C = linspace(0.7, 1.2, 4);  %Range of C values to test
pulsewidth = linspace(10, 70, 5); %Range of pw values to test
A = linspace(0.1, 1, 3); %Range of A values being tested
lambda = linspace(725, 755, 5)*2*10^-9;  %Range of lambda values to test

figure

subplot(3, 4, 1)
pcolor(lambda, pulsewidth, squeeze(Purity(1, :, 1, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 0.1 kW, C = 70 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 2)
pcolor(lambda, pulsewidth, squeeze(Purity(2, :, 1, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 0.1 kW, C = 86.7 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 3)
pcolor(lambda, pulsewidth, squeeze(Purity(3, :, 1, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 0.1 kW, C = 103.3 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 4)
pcolor(lambda, pulsewidth, squeeze(Purity(4, :, 1, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 0.1 kW, C = 120 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 5)
pcolor(lambda, pulsewidth, squeeze(Purity(1, :, 2, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 0.55 kW, C = 70 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 6)
pcolor(lambda, pulsewidth, squeeze(Purity(2, :, 2, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 0.55 kW, C = 86.7 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 7)
pcolor(lambda, pulsewidth, squeeze(Purity(3, :, 2, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 0.55 kW, C = 103.3 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 8)
pcolor(lambda, pulsewidth, squeeze(Purity(4, :, 2, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 0.55 kW, C = 120 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 9)
pcolor(lambda, pulsewidth, squeeze(Purity(1, :, 3, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 1 kW, C = 70 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 10)
pcolor(lambda, pulsewidth, squeeze(Purity(2, :, 3, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 1 kW, C = 86.7 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 11)
pcolor(lambda, pulsewidth, squeeze(Purity(3, :, 3, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title(['A = 1 kW, C = 103.3 1/m'], 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

subplot(3, 4, 12)
pcolor(lambda, pulsewidth, squeeze(Purity(4, :, 3, :)))
shading interp
xlabel('\lambda_F / m', 'fontsize', 15)
ylabel('\tau / ps', 'fontsize', 15)
title('A = 1 kW, C = 120 1/m', 'fontsize', 15)
clim([0 0.3785])
set(gca, 'FontSize', 12);

colorbar('Position', [0.92, 0.1, 0.02, 0.8]);