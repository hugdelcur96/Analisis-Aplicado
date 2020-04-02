clear;   close all;   clc;


% level sets
f=@(X,Y) 0.26 * (X.^2 + Y.^2) - 0.48 .*X.*Y;
%x0 = [1.5;1];
x0=[-1.5;-1];
delta = 1;
stepsize =  0.01;  % mas chico hace los conjuntos de nivel mas detallados
[X,Y] = meshgrid(-2:stepsize:2);
z = f(X,Y);
niveles = [-2:0.1:2];
contour(X,Y,z, niveles)

axis equal

f = @(x) (0.26 * (x(1)^2 + x(2)^2) - 0.48 * x(1) * x(2));;
gk = apGrad(f,x0);
Bk = apHess(f,x0);
pN = -inv(Bk)*gk;
pC = pCauchy(Bk,gk,delta);
pDog = pDogLeg(Bk,gk,delta);

hold on
dirN = quiver(x0(1), x0(2),pN(1), pN(2),'Color', 'green',   'LineWidth', 2);
dirC = quiver(x0(1), x0(2), pC(1), pC(2),'Color', 'red',   'LineWidth', 2);
dirDog = quiver(x0(1), x0(2), pDog(1), pDog(2),'Color', 'blue', 'LineWidth', 2);
menosGrad = quiver(x0(1), x0(2), gk(1), gk(2),'Color', 'black', 'LineWidth', 2);
regionConfianza = viscircles(x0', delta,'LineStyle', '--', 'Color', 'm', 'LineWidth',1);
legend([dirN, dirC, dirDog, regionConfianza], {'Dirección Newton', ... 
    'Dirección Cauchy', 'Dirección Dogleg', 'Región de confianza'});
hold off
grid on
