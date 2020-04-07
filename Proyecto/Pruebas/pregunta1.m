clear;   close all;   clc;

%% Define f with two arguments  and  so that
% it can be evaluated for matrices of values X, Y
f = @(X, Y)  (0.26) * (X.^2 + Y.^2) - 0.48 .* X .* Y;
Df = @(X, Y) [0.52 .* X - 0.48 .* Y; 0.52 .* Y - 0.48 .* X];
D2f = [0.52 -0.48; -0.48 0.52];


%% define point and trust region radius
x0    = [0.25;0.25];
delta = 0.5;


%% plot f in cartesian coordinates arround x0
showPlot = false;
if showPlot
	Delta = 1.1*delta;
	uniGrid = linspace(x0(1)-Delta, x0(1)+Delta, 32);
	[X,Y]   = meshgrid(uniGrid, uniGrid);
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
m  = @(X, Y)   f(X, Y);
% evaluation and plot
Z  = m(X,Y);
s1 = mesh(X,Y,Z);
view(130, 15)
hold off