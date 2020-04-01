clear;   close all;   clc;

f = @(x) (0.26 * (x(1)^2 + x(2)^2) - 0.48 * x(1) * x(2));
x0 = [0.5; 0.5];
itmax = 100;

[xk, msg, iter] = mRC1(f, x0, itmax)