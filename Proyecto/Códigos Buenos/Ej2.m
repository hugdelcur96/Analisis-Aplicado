clear;   close all;   clc;

% level sets
% f=@(X,Y) 0.26 * (X.^2 + Y.^2) - 0.48 .*X.*Y;
f=@(X,Y) 0.125 - 0.3.*X + 0.2.*Y + 0.26 .* X.^2 - 0.48 .* X .* Y + 0.26 .* Y.^2; %Modelo cuadratico en x0
%x0 = [1.5;1];
x0=[-1.5;-1];
delta = 1;
stepsize =  0.01;  % mas chico hace los conjuntos de nivel mas detallados
[X,Y] = meshgrid(-2.5:stepsize:2.5);
z = f(X,Y);
niveles = [-5:0.1:5];
contour(X,Y,z, niveles)

% Para el circulo
theta = linspace(0, 2*pi, 100);
xcirc = delta * cos(theta) + x0(1);
ycirc = delta * sin(theta) + x0(2);


axis equal

f1=fmatyas();
gk = apGrad(f1,x0);
Bk = apHess(f1,x0);
pN = -inv(Bk)*gk;
pC = pCauchy(Bk,gk,delta);
pDog = pDogLeg(Bk,gk,delta);

hold on
dirN = quiver(x0(1), x0(2),pN(1), pN(2),'Color', 'green',   'LineWidth', 2);
dirC = quiver(x0(1), x0(2), pC(1), pC(2),'Color', 'red',   'LineWidth', 2);
dirDog = quiver(x0(1), x0(2), pDog(1), pDog(2),'Color', 'blue', 'LineWidth', 2);
regionConfianza = plot(xcirc, ycirc, 'r--');
%regionConfianza = viscircles(x0', delta,'LineStyle', '--', 'Color', 'm', 'LineWidth',1);
legend([dirN, dirC, dirDog, regionConfianza], {'Dirección Newton', ... 
    'Dirección Cauchy', 'Dirección Dogleg', 'Región de confianza'});
hold off
grid on
