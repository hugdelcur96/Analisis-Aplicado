clear;   close all;   clc;

[f, Df, D2f] = fpascal();
x0 = 4*ones(4,1);
itmax = 1000;


[xk, msg] = mRC1(f, x0, itmax)