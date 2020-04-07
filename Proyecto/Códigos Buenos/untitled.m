clear all
clc

f=@(x)((1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*(30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2)));
x0 = [0.5; -0.5];
itmax = 1000;

[x, msg, iter] = mRC1(f, x0, itmax)
[x, msg, iter] = mRC2(f, x0, itmax)