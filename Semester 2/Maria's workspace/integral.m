clear all, close all
clc

% parameters - placeholders atm
B = 1; 
p = 3; 

% grid domain vals
z = -3:.1:3; 

% function
f = p * exp(1i * B * z);

I = trapz(f);
disp(I);
