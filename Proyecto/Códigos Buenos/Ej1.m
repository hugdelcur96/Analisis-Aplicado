clear;   close all;   clc;

f = @(X,Y) 0.26 * (X.^2 + Y.^2) - 0.48 .*X.*Y;

%% define point and trust region radius
x0    = [0.5;0.5];
delta = 1;

%% plot f in cartesian coordinates arround x0
showPlot = true;
if showPlot
	Delta = 1.1*delta;
	uniGrid = linspace(x0(1)-Delta, x0(1)+Delta, 32);
	[X,Y] = meshgrid(uniGrid, uniGrid);
	Z  = f(X,Y);
	s2 = surf(X,Y,Z);
end


%% plot quadratic model in trust region with polar coordinates arround x0
hold on
% polar coordinates arround x0
[T,R] = meshgrid(linspace(0,2*pi,64),linspace(0,delta,16));
X     = R.*cos(T) +x0(1);
Y     = R.*sin(T) +x0(2);
% the quadratic model (simple in this case)
m  = @(X,Y) 0.01 - 0.02.*X + 0.02.*Y +  0.26.*X.^2 - 0.48 .* X .* Y + 0.26 .* Y.^2;
% evaluation and plot
Z  = m(X,Y);
s1 = mesh(X,Y,Z);
view(130, 15)
hold off
