clear;   close all;   clc;

% level sets
f=@(X,Y) 0.26 * (X.^2 + Y.^2) - 0.48 .*X.*Y;
%x0 = [1.5;1];
x0=[3;4];
delta = 1;
stepsize =  0.01;  % mas chico hace los conjuntos de nivel mas detallados
[X,Y] = meshgrid(-2:stepsize:5);
z = f(X,Y);
niveles = [-2:0.1:5];
contour(X,Y,z, niveles)

axis equal
f=fmatyas();
it = 1000;
x=zeros(129,1);
y=zeros(129,1);
x(1)=x0(1);
y(1)=x0(2);
j=2;
for i=1:127
    [x1, msg1, h] = mRC1(f, x0, i);
    x(j)=x1(1);
    y(j)=x1(2);
    j=j+1;
end
hold on
x = x;
y = y;
plot(x,y,'--d')
hold off

