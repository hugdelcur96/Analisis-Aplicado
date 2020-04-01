clear all
clc

f = @(x) [0.5, 1] * x'.^2;
Df = @(x) [x(1), 2 * x(2)]';
D2f = [1 0; 0 2];

x0 = [1, 2];
g = Df(x0);
B = D2f;
delta = 2;

p = pDogLeg(B, g, delta);