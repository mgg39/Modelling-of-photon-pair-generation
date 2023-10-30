clear all
clc

syms F z t r2 S k s1 s2 g P p1 p2

UF = i*(diff(F, z)) == -r2*(diff(F, t, t)) - S*conj(F)*exp(i*k*z);
US = i*(diff(S, z))== -s1*(diff(S, t)) - s2*(diff(S, t, t)) - ((F^2)/2)*exp(-i*k*z) +g*P;
UP = i*(diff(P, z))== -p1*(diff(P, t)) - p2*(diff(P, t, t)) + g*S;