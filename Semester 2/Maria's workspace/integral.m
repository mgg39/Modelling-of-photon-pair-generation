clear all, close all
clc

% parameters - placeholders atm
B = 1; 
O = 2; 
p = 3; 

% grid domain vals
z = -3:.1:3; 

% function
f = p * exp(1i * B * O * z);

I = trapz(f);
disp(I);
