%La función Matyas es convexa
%Con mínimo en (0,0)
close all; clc;
format long

%deltaMax = 1.5;
%Definimos x0 de tal manera que la distancia con xmin = (0,0) mayor q
%3*deltaMax=4.5
xmin = [0;0];
x0 = [3;4];
f =  fmatyas();
%Definimos el modelo cuadrático: necesitamos la hessiana Bk y el gradiente
%g
gk = apGrad(f,x0);
Bk = apHess(f,x0);

% Usaremos un numero de iteraciones muy grande para que no influya en ver que metodo se
% aproxima mas rapido.

it = 1000;

[x1, msg1, iter1] = mRC1(f, x0, it);
[x2, msg2, iter2] = mRC2(f, x0, it);
%Veamos cuantas iteraciones se realizo para cada proceso
fprintf('\n\t%s \t\t\t%s\n' ,'Método', 'Iteraciones');
fprintf('\t--------------------------------------------');

fprintf('\n\t%s  \t\t%d', 'pCauchy', iter1);
fprintf('\n\t%s  \t\t%d', 'pDogleg', iter2);
fprintf('\n');
fprintf('\t--------------------------------------------');
fprintf('\n');

%Ahora veremos los errores, la norma de la derivada y el punto para las
%ultimas 8 iteraciones

y=zeros(8,2);
z=zeros(8,2);
der1=zeros(8,1);
der2=zeros(8,1);
fx1=zeros(8,1);
fx2=zeros(8,1);
err1=zeros(8,1);
err2=zeros(8,1);
j=1;
for i=iter1-8:iter1
    [x, msg1, h] = mRC1(f, x0, i);
    y(j,1)=x(1);
    y(j,2)=x(2);
    der1(j)=norm(apGrad(f,x));
    fx1(j)=(0.26 * (x(1)^2 + x(2)^2) - 0.48 * x(1) * x(2));
    err1(j) = norm(x-xmin);
    j=j+1;
end
j=1;
for i=iter2-4:iter2
    [x, msg2, m] = mRC2(f, x0, i);
    z(j,1)=x(1);
    z(j,2)=x(2);
    der2(j)=norm(apGrad(f,x));
    fx2(j)=(0.26 * (x(1)^2 + x(2)^2) - 0.48 * x(1) * x(2));
    err2(j) = norm(x-xmin);
    j=j+1;
end

%Genera la Tabla para mRC1
fprintf('Tabla para mRC1')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t%s \t\t%s \t\t\t%s\n' ,'Iteracion','Valor X','Valor Y','NormGrad','Evaluada','Error');
fprintf('\t-------------------------------------------------------------------------------------------------');

fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '107',y(2,1),y(2,2),der1(2),fx1(2),err1(2));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '108',y(3,1),y(3,2),der1(3),fx1(3),err1(3));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '109',y(4,1),y(4,2),der1(4),fx1(4),err1(4));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '110',y(5,1),y(5,2),der1(5),fx1(5),err1(5));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '111',y(6,1),y(6,2),der1(6),fx1(6),err1(6));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '112',y(7,1),y(7,2),der1(7),fx1(7),err1(7));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '113',y(8,1),y(8,2),der1(8),fx1(8),err1(8));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '114',y(9,1),y(9,2),der1(9),fx1(9),err1(9));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');

%Genera la Tabla para mRC2
fprintf('Tabla para mRC2')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t%s \t\t%s \t\t\t%s\n' ,'Iteracion','Valor X','Valor Y','NormGrad','Evaluada','Error');
fprintf('\t-------------------------------------------------------------------------------------------------');

fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '1',z(2,1),z(2,2),der2(1),fx2(2),err2(2));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '2',z(3,1),z(3,2),der2(2),fx2(3),err2(3));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '3',z(4,1),z(4,2),der2(3),fx2(4),err2(4));
fprintf('\n\t%s \t\t%d \t\t%d \t\t%d \t%d \t\t%.4d', '4',z(5,1),z(5,2),der2(5),fx2(5),err2(5));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');