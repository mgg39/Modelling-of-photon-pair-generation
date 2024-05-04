clear all, close all
clc

%% 

N = 2048;

x = linspace(-5, 5, N);
y = linspace(-30, 30, N);

gaus = exp(-(x.^2));
Sinc = sinc(y, 2);

[G, S] = meshgrid(gaus, Sinc);

pattern = G.*S.*1;

figure
subplot(4, 3, [2,3,5,6,8,9])
pcolor(x, y, pattern)
shading interp
fontsize(gca, 18, 'points')
title("F(x,y) = X(x)Y(y)", "fontsize", 20)
subplot(4, 3, [1,4,7])
plot(Sinc,y)
xlabel("sinc(y)", "fontsize", 20)
ylabel("y", "rotation", 0, "fontsize", 20)
fontsize(gca, 18, 'points')
title("Y(y)", "fontsize", 20)
subplot(4, 3, [11,12])
plot(x, gaus)
xlabel("x", "fontsize", 20)
ylabel("exp(-(x^2))", "rotation", 0, "position", [-6.3 0.4], "fontsize", 20)
fontsize(gca, 18, 'points')
title("X(x)", "fontsize", 20)
