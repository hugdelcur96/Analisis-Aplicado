clear all
clc

a = 1
b = 1
c = 2

f = @(x) (a*x(1)^2 + b*x(2)^2 + c*x(3)^2) 

x0 = [1; 1; 1]

B = 2*[a 0 0; 0 b 0; 0 0 c]
g = 2*[a * x0(1); b * x0(2); c * x0(3)]
pB = - inv(B) * g
posdef = dot(g, B * g)
normapB = norm(pB)
normag = norm(g)
pU = -normag^2 / posdef * g

normag^3 / posdef < normapB

liminf = normag^3 / posdef
limsup = normapB

delta =1.5

a1 = norm(pB - pU)^2
b1 = 2 * dot(pU, pB - pU)
c1 = norm(pU)^2 - delta^2
alfa1 = (-b1 + sqrt(b1^2 - 4*a1*c1)) / (2 * a1)
alfa2 = (-b1 - sqrt(b1^2 - 4*a1*c1)) / (2 * a1)

pDL = pU + alfa1 * (pB - pU)

pDLcomunidad = pDogLeg(B, g, delta)
            
            
            
            
            